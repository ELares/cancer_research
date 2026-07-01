# Reproducibility container for ELares/cancer_research (#542, the unshipped half of #350).
#
# Pins Python + the Rust toolchain so the project's checked results reproduce
# end-to-end with one command. Build then run:
#
#     docker build -t cancer-research .
#     docker run --rm cancer-research          # == `make reproduce`
#
# The default command runs `make reproduce` (pytest + the Rust workspace test
# suite, which includes the calibration-regression gate that re-checks every
# data-anchored leg, plus a tangible headline regeneration).
FROM python:3.12-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential curl ca-certificates git pkg-config make \
    && rm -rf /var/lib/apt/lists/*

# Pinned Rust toolchain — matches simulations/rust-toolchain.toml (1.96.0).
ENV RUSTUP_HOME=/usr/local/rustup \
    CARGO_HOME=/usr/local/cargo \
    PATH=/usr/local/cargo/bin:$PATH
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | \
        sh -s -- -y --no-modify-path --default-toolchain 1.96.0 --profile minimal

WORKDIR /work
COPY . .

# Pinned environment tooling. uv installs the committed lockfile as-is; it does
# not re-resolve here. maturin stays available for optional local bindings work.
RUN pip install --no-cache-dir uv==0.11.26 maturin==1.14.1 && \
    uv sync --locked --group test

ENV PATH=/work/.venv/bin:$PATH

# Pre-build the Rust workspace so `make reproduce` does not pay the cold build.
RUN cd simulations && cargo build --release --workspace

CMD ["make", "reproduce"]
