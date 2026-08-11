#!/usr/bin/env bash
# Launch a one-shot EC2 box for the 1e11 tier of the rare-event sweep, with
# teardown that does not depend on anyone remembering.
#
# WHY THIS EXISTS AT ALL
# ----------------------
# Only the 1e11 tier needs a cloud box. Measured on a 14-core laptop, 1e9 is 40
# minutes and 1e10 is an overnight run; 1e11 across three conditions is 639
# CPU-hours, which is 70 hours locally and about 10 on 64 vCPUs. Everything
# smaller is cheaper to run at home than to orchestrate.
#
# FOUR INDEPENDENT KILL SWITCHES, because the failure that costs money is the
# one where the job dies and the instance does not:
#
#   1. --instance-initiated-shutdown-behavior terminate
#      so any `shutdown` from inside terminates rather than stopping. A stopped
#      instance still bills for its EBS.
#   2. a dead-man `shutdown -h +$MAX_MINUTES` scheduled as the FIRST action of
#      user-data, before the job starts. If the job hangs, wedges, or the build
#      fails, the box still dies on schedule.
#   3. the job script calls `shutdown -h now` when it finishes, success or
#      failure, so the common case does not wait for the dead-man.
#   4. a CloudWatch alarm terminating on sustained low CPU, which catches the
#      case where the job exits without running its trap and the box idles.
#
# WHY THE ROOT VOLUME OUTLIVES THE INSTANCE
# ------------------------------------------
# This account has no S3 access, so there is nowhere durable to stream results.
# The root volume is therefore created with DeleteOnTermination=false: when the
# dead-man fires, the results survive on the volume and can be recovered by
# attaching it elsewhere. Without that, an auto-terminate would destroy exactly
# the thing the run existed to produce. The volume is tagged so it is findable,
# and deleting it afterwards is a manual step -- deliberately, because an
# automatic delete would reintroduce the problem.
#
# The sweep is also incremental: every (condition, scale) result is appended as
# soon as it completes and is never recomputed, so an interrupted run loses at
# most one cell.
#
# PRIVACY: nothing account-specific is hardcoded here. Region, key, subnet and
# AMI come from the environment or from a lookup at run time, so this file is
# safe in a public repository.
#
# Usage:
#   AWS_REGION=us-west-2 KEY_NAME=my-key ./scripts/aws/launch_scale_run.sh
#   DRY_RUN=1 ... to print the plan without launching

set -euo pipefail

REGION="${AWS_REGION:-$(aws configure get region)}"
INSTANCE_TYPE="${INSTANCE_TYPE:-c7g.16xlarge}"   # 64 vCPU Graviton
MAX_MINUTES="${MAX_MINUTES:-900}"                # 15 h dead-man; job needs ~10
REPO_URL="${REPO_URL:-https://github.com/ELares/cancer_research.git}"
REPO_REF="${REPO_REF:-main}"
TAG_NAME="${TAG_NAME:-ferro-scale-1e11}"
DRY_RUN="${DRY_RUN:-0}"

if [[ -z "${KEY_NAME:-}" ]]; then
  echo "set KEY_NAME to an EC2 key pair you hold the private half of" >&2
  exit 2
fi

# Latest Amazon Linux 2023 for arm64, resolved rather than pinned so this does
# not rot; pin it if a run must be reproducible to the AMI.
AMI="${AMI:-$(aws ssm get-parameters --region "$REGION" \
  --names /aws/service/ami-amazon-linux-latest/al2023-ami-kernel-default-arm64 \
  --query 'Parameters[0].Value' --output text)}"

SUBNET="${SUBNET_ID:-$(aws ec2 describe-subnets --region "$REGION" \
  --filters Name=default-for-az,Values=true \
  --query 'Subnets[0].SubnetId' --output text)}"

USER_DATA=$(cat <<EOF
#!/bin/bash
# (2) DEAD-MAN FIRST. Before the build, before the job, before anything that can
# fail. Every later step is optional; this one is not.
shutdown -h +${MAX_MINUTES}

exec > >(tee -a /var/log/ferro-run.log) 2>&1
set -x

dnf -y install git gcc tmux
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
source /root/.cargo/env

cd /root
git clone --depth 1 --branch ${REPO_REF} ${REPO_URL} repo
cd repo/simulations
cargo build --release -p sim-scale

# The gate: an engine that cannot reproduce the committed figure does not get
# to spend ten hours. If this fails we shut down immediately rather than
# producing numbers nobody should trust.
if ! ./target/release/sim-scale --verify; then
  echo "SELF-CHECK FAILED — aborting before the expensive run"
  shutdown -h now
  # EXIT, do not merely schedule the halt. \`shutdown\` is asynchronous: it
  # returns immediately and the system takes time to come down, so without this
  # the script falls straight through to the expensive sweep and runs it during
  # the halt window -- having just printed "aborting". Observed for real: an
  # accidental local execution of this body hit a failing self-check, printed
  # the abort message, had \`shutdown\` refuse with "NOT super-user", and then
  # started a 1e11 sweep anyway. The gate did not gate.
  exit 1
fi

cd /root/repo
mkdir -p /root/results
OUT=/root/results/rare-event-sweep-1e11.jsonl

# TWO PASSES, deliberately. A condition that dies to a signal is recorded and
# skipped rather than retried within a pass, and the sweep is resumable by
# design: a completed (condition, n) is never recomputed, so a second
# invocation costs nothing for the ones that worked and retries only the ones
# that did not. A local 1e10 run lost two good conditions to one dead one
# before the driver was fixed to continue; here a lost condition would also
# mean a whole second instance, so it is worth one free retry.
for pass in 1 2; do
  echo "=== sweep pass \$pass ==="
  python3 scripts/rare_event_sweep.py --scales 1e11 --out "\$OUT" && break
  echo "pass \$pass left conditions incomplete; retrying the remainder"
done

# Say plainly what came back, so the console log alone answers "did this work"
# without attaching the volume.
echo "=== results ==="
wc -l < "\$OUT" 2>/dev/null || echo "NO RESULTS FILE"
cat "\$OUT" 2>/dev/null

# (3) done, so stop paying immediately rather than waiting for the dead-man
shutdown -h now
EOF
)

echo "plan:"
echo "  region        $REGION"
echo "  instance      $INSTANCE_TYPE"
echo "  dead-man      shutdown -h +${MAX_MINUTES} min, set before the job starts"
echo "  on shutdown   TERMINATE (not stop)"
echo "  root volume   DeleteOnTermination=false, so results survive the teardown"
echo "  tag           Name=$TAG_NAME"
if [[ "$DRY_RUN" == "1" ]]; then
  echo "DRY_RUN=1, not launching"
  exit 0
fi

IID=$(aws ec2 run-instances --region "$REGION" \
  --image-id "$AMI" --instance-type "$INSTANCE_TYPE" \
  --key-name "$KEY_NAME" --subnet-id "$SUBNET" \
  --instance-initiated-shutdown-behavior terminate \
  --block-device-mappings "[{\"DeviceName\":\"/dev/xvda\",\"Ebs\":{\"VolumeSize\":40,\"DeleteOnTermination\":false}}]" \
  --tag-specifications \
    "ResourceType=instance,Tags=[{Key=Name,Value=$TAG_NAME},{Key=Purpose,Value=rare-event-sweep},{Key=AutoTerminate,Value=true}]" \
    "ResourceType=volume,Tags=[{Key=Name,Value=$TAG_NAME-results}]" \
  --user-data "$USER_DATA" \
  --query 'Instances[0].InstanceId' --output text)

echo "launched $IID"

# (4) last resort: terminate if the box sits idle, which means the job is gone
aws cloudwatch put-metric-alarm --region "$REGION" \
  --alarm-name "${TAG_NAME}-idle-terminate" \
  --alarm-description "terminate ${TAG_NAME} if CPU stays low, i.e. the job has ended" \
  --metric-name CPUUtilization --namespace AWS/EC2 --statistic Average \
  --period 300 --evaluation-periods 6 --threshold 5 \
  --comparison-operator LessThanThreshold \
  --dimensions "Name=InstanceId,Value=$IID" \
  --alarm-actions "arn:aws:automate:${REGION}:ec2:terminate" >/dev/null

echo "idle-terminate alarm armed (30 min under 5% CPU)"
echo
echo "watch:      aws ec2 describe-instances --region $REGION --instance-ids $IID --query 'Reservations[0].Instances[0].State.Name' --output text"
echo "results on: the tagged volume, which survives termination"
echo "kill now:   aws ec2 terminate-instances --region $REGION --instance-ids $IID"
