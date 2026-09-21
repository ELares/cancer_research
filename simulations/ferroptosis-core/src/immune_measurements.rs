//! Passive observations of the canonical immune comparison.
//!
//! This module has no RNG and never receives mutable simulation state. A cell
//! identity is a lattice index only because the measured arms forbid regrowth
//! and additional killing routes. Every floating-point aggregate uses stable
//! cell-index order, independently of the simulation's Rayon scheduling.

use crate::grid::GridCell;
use serde::Serialize;

#[derive(Clone, Debug, Serialize)]
pub struct FerroptoticEvent {
    cell_index: usize,
    death_step: u32,
    death_lp: f64,
    scheduled_release_step: u32,
    release_step: Option<u32>,
    release_lp: Option<f64>,
    release_damp: Option<f64>,
    horizon_lp: Option<f64>,
    terminal_damp: Option<f64>,
}

#[derive(Clone, Debug, Serialize)]
pub struct ImmuneKillEvent {
    cell_index: usize,
    step: u32,
    local_damp: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct EligibleCell {
    cell_index: usize,
    first_step: u32,
    last_step: u32,
    opportunities: usize,
    local_damp_sum: f64,
    #[serde(flatten, skip_serializing_if = "Option::is_none")]
    activation: Option<Activation>,
}

#[derive(Clone, Debug, Serialize)]
pub struct Step {
    step: u32,
    ferroptotic_deaths: usize,
    completed_releases: usize,
    released_damp: f64,
    eligible_cells: usize,
    eligible_local_damp_sum: f64,
    immune_kills: usize,
    #[serde(flatten, skip_serializing_if = "Option::is_none")]
    activation: Option<Activation>,
}

/// Opportunity-weighted activation diagnostics, sampled before each immune draw.
/// These optional fields are omitted entirely from the historical 3D v1 ledger.
#[derive(Clone, Debug, Default, Serialize)]
pub struct Activation {
    activation_sum: f64,
    max_local_damp: Option<f64>,
    damp_ge_kd_opportunities: usize,
    damp_ge_9kd_opportunities: usize,
}

impl Activation {
    fn observe(&mut self, damp: f64, kd: f64) {
        self.activation_sum += damp / (damp + kd);
        self.max_local_damp = Some(
            self.max_local_damp
                .map_or(damp, |previous| previous.max(damp)),
        );
        self.damp_ge_kd_opportunities += usize::from(damp >= kd);
        self.damp_ge_9kd_opportunities += usize::from(damp >= 9.0 * kd);
    }
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct Terminal {
    censored_deaths: usize,
    terminal_additions: usize,
    terminal_damp: f64,
    damp_before_terminal: f64,
    damp_after_terminal: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct Measurements {
    pub ferroptotic_events: Vec<FerroptoticEvent>,
    pub immune_kill_events: Vec<ImmuneKillEvent>,
    pub eligible_cells: Vec<EligibleCell>,
    pub steps: Vec<Step>,
    pub terminal: Terminal,
    #[serde(skip)]
    event_by_cell: Vec<Option<usize>>,
    #[serde(skip)]
    eligibility_by_cell: Vec<Option<EligibleCell>>,
    #[serde(skip)]
    eligible_this_step: Vec<usize>,
    #[serde(skip)]
    post_death_steps: u32,
    #[serde(skip)]
    damp_per_lp: f64,
    #[serde(skip)]
    immune_start_step: u32,
    #[serde(skip)]
    damp_kill_threshold: f64,
    #[serde(skip)]
    activation_kd: Option<f64>,
}

/// Match the production admission checks, including equality at the threshold.
/// The delay lives here as well because a DAMP-positive cell before that delay
/// has no immune-draw opportunity. Preserve the original `<` comparison rather
/// than changing its NaN semantics to `>=`.
fn eligible(
    is_tumor: bool,
    dead: bool,
    local_damp: f64,
    step: u32,
    start: u32,
    threshold: f64,
) -> bool {
    if step < start || dead || !is_tumor {
        return false;
    }
    if local_damp < threshold {
        return false;
    }
    true
}

impl Measurements {
    pub fn new(
        n_cells: usize,
        post_death_steps: u32,
        damp_per_lp: f64,
        immune_start_step: u32,
        damp_kill_threshold: f64,
        activation_kd: Option<f64>,
    ) -> Self {
        assert!(activation_kd.is_none_or(|kd| kd.is_finite() && kd > 0.0));
        assert!(
            post_death_steps > 0,
            "the canonical release delay is positive"
        );
        Self {
            ferroptotic_events: Vec::new(),
            immune_kill_events: Vec::new(),
            eligible_cells: Vec::new(),
            steps: Vec::new(),
            terminal: Terminal::default(),
            event_by_cell: vec![None; n_cells],
            eligibility_by_cell: vec![None; n_cells],
            eligible_this_step: Vec::new(),
            post_death_steps,
            damp_per_lp,
            immune_start_step,
            damp_kill_threshold,
            activation_kd,
        }
    }

    /// Called immediately after the biochemical loop, before other death
    /// routes or iron diffusion. Death LP is still its true threshold-crossing
    /// value; the stored grace-end value is the value actually injected this
    /// step, before diffusion and clearance.
    pub fn after_biochemistry(&mut self, cells: &[GridCell], step: u32) {
        assert_eq!(self.steps.len(), step as usize);
        let mut row = Step {
            step,
            ferroptotic_deaths: 0,
            completed_releases: 0,
            released_damp: 0.0,
            eligible_cells: 0,
            eligible_local_damp_sum: 0.0,
            immune_kills: 0,
            activation: self.activation_kd.map(|_| Activation::default()),
        };
        for (idx, gc) in cells.iter().enumerate() {
            if !gc.is_tumor || !gc.state.dead {
                continue;
            }
            let Some(ds) = gc.state.death_step else {
                continue; // Immune deaths deliberately have no death_step.
            };
            if ds == step {
                assert!(
                    self.event_by_cell[idx].is_none(),
                    "cell identity was reused"
                );
                self.event_by_cell[idx] = Some(self.ferroptotic_events.len());
                self.ferroptotic_events.push(FerroptoticEvent {
                    cell_index: idx,
                    death_step: ds,
                    death_lp: gc.state.lp,
                    scheduled_release_step: ds + self.post_death_steps,
                    release_step: None,
                    release_lp: None,
                    release_damp: None,
                    horizon_lp: None,
                    terminal_damp: None,
                });
                row.ferroptotic_deaths += 1;
            }
            if step == ds + self.post_death_steps {
                let event = &mut self.ferroptotic_events
                    [self.event_by_cell[idx].expect("release must have an observed death")];
                assert!(event.release_step.is_none(), "duplicate release");
                let damp = gc.lp_at_grace_end * self.damp_per_lp;
                event.release_step = Some(step);
                event.release_lp = Some(gc.lp_at_grace_end);
                event.release_damp = Some(damp);
                row.completed_releases += 1;
                row.released_damp += damp;
            }
        }
        self.steps.push(row);
    }

    /// Called after DAMP diffusion and immediately before immune killing.
    pub fn before_immune(&mut self, cells: &[GridCell], damp: &[f64], step: u32) {
        assert_eq!(cells.len(), damp.len());
        self.eligible_this_step.clear();
        let row = self.steps.last_mut().expect("biochemistry observed first");
        assert_eq!(row.step, step);
        for (idx, gc) in cells.iter().enumerate() {
            if !eligible(
                gc.is_tumor,
                gc.state.dead,
                damp[idx],
                step,
                self.immune_start_step,
                self.damp_kill_threshold,
            ) {
                continue;
            }
            self.eligible_this_step.push(idx);
            row.eligible_cells += 1;
            row.eligible_local_damp_sum += damp[idx];
            let cell = self.eligibility_by_cell[idx].get_or_insert(EligibleCell {
                cell_index: idx,
                first_step: step,
                last_step: step,
                opportunities: 0,
                local_damp_sum: 0.0,
                activation: self.activation_kd.map(|_| Activation::default()),
            });
            cell.last_step = step;
            cell.opportunities += 1;
            cell.local_damp_sum += damp[idx];
            if let Some(kd) = self.activation_kd {
                cell.activation
                    .as_mut()
                    .expect("activation enabled")
                    .observe(damp[idx], kd);
                row.activation
                    .as_mut()
                    .expect("activation enabled")
                    .observe(damp[idx], kd);
            }
        }
    }

    /// Only the immune loop ran since before_immune; an eligible cell that is
    /// now dead therefore corresponds to one actual successful immune draw.
    pub fn after_immune(&mut self, cells: &[GridCell], damp: &[f64], step: u32) {
        let row = self.steps.last_mut().expect("biochemistry observed first");
        assert_eq!(row.step, step);
        for &idx in &self.eligible_this_step {
            if cells[idx].state.dead {
                self.immune_kill_events.push(ImmuneKillEvent {
                    cell_index: idx,
                    step,
                    local_damp: damp[idx],
                });
                row.immune_kills += 1;
            }
        }
    }

    /// The horizon flush is an addition AFTER all immune updates. It is
    /// recorded separately from completed releases, even when its intended
    /// release step equals the horizon exactly.
    pub fn before_terminal(&mut self, cells: &[GridCell], damp: &[f64], n_steps: u32) {
        assert_eq!(self.steps.len(), n_steps as usize);
        self.terminal.damp_before_terminal = damp.iter().sum();
        for event in &mut self.ferroptotic_events {
            if event.scheduled_release_step >= n_steps {
                assert!(event.release_step.is_none());
                let gc = &cells[event.cell_index];
                assert!(gc.is_tumor && gc.state.dead);
                let addition = gc.state.lp * self.damp_per_lp;
                event.horizon_lp = Some(gc.state.lp);
                event.terminal_damp = Some(addition);
                self.terminal.censored_deaths += 1;
                self.terminal.terminal_additions += 1;
            } else {
                assert!(event.release_step.is_some(), "missing completed release");
            }
        }
        // Use index order, matching the actual terminal flush, rather than
        // event order (death step, then index).
        self.terminal.terminal_damp = self
            .event_by_cell
            .iter()
            .flatten()
            .filter_map(|&i| self.ferroptotic_events[i].terminal_damp)
            .sum();
        self.eligible_cells = self.eligibility_by_cell.iter().flatten().cloned().collect();
    }

    pub fn after_terminal(&mut self, damp: &[f64]) {
        self.terminal.damp_after_terminal = damp.iter().sum();
    }

    pub fn validate_totals(
        &self,
        ferroptosis_kills: usize,
        immune_kills: usize,
        total_dead: usize,
        total_damp: f64,
    ) {
        let ferro = self.ferroptotic_events.len();
        let immune = self.immune_kill_events.len();
        assert_eq!(ferro, ferroptosis_kills);
        assert_eq!(immune, immune_kills);
        assert_eq!(ferro + immune, total_dead);
        assert_eq!(
            ferro,
            self.steps
                .iter()
                .map(|s| s.completed_releases)
                .sum::<usize>()
                + self.terminal.censored_deaths
        );
        assert!(immune <= self.eligible_cells.len());
        assert_eq!(
            self.eligible_cells
                .iter()
                .map(|c| c.opportunities)
                .sum::<usize>(),
            self.steps.iter().map(|s| s.eligible_cells).sum::<usize>()
        );
        assert_eq!(self.terminal.damp_after_terminal, total_damp);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grid::TumorGrid3D;
    const IMMUNE_START_STEP: u32 = 60;
    const DAMP_KILL_THRESHOLD: f64 = 0.01;

    fn eligible(is_tumor: bool, dead: bool, local_damp: f64, step: u32) -> bool {
        super::eligible(
            is_tumor,
            dead,
            local_damp,
            step,
            IMMUNE_START_STEP,
            DAMP_KILL_THRESHOLD,
        )
    }

    #[test]
    fn eligibility_uses_delay_liveness_and_exact_threshold() {
        assert!(!eligible(true, false, 100.0, IMMUNE_START_STEP - 1));
        assert!(eligible(
            true,
            false,
            DAMP_KILL_THRESHOLD,
            IMMUNE_START_STEP
        ));
        assert!(!eligible(
            true,
            false,
            DAMP_KILL_THRESHOLD / 2.0,
            IMMUNE_START_STEP
        ));
        assert!(!eligible(false, false, 100.0, IMMUNE_START_STEP));
        assert!(!eligible(true, true, 100.0, IMMUNE_START_STEP));
    }

    #[test]
    fn release_boundary_and_horizon_censoring_are_distinct() {
        let mut grid = TumorGrid3D::generate(3, 3, 3, 20.0, 42);
        let cells = &mut grid.cells;
        // Artificial events isolate the observer's timing from stochastic
        // biochemistry. One releases at the last simulated step, the other
        // exactly at the horizon and must remain censored.
        for cell in cells.iter_mut() {
            cell.is_tumor = false;
        }
        let mut observer = Measurements::new(
            cells.len(),
            2,
            1.5,
            IMMUNE_START_STEP,
            DAMP_KILL_THRESHOLD,
            None,
        );
        for step in 0..3 {
            if step < 2 {
                let cell = &mut cells[step as usize];
                cell.is_tumor = true;
                cell.state.dead = true;
                cell.state.death_step = Some(step);
                cell.state.lp = 11.0 + step as f64;
            }
            if step == 2 {
                cells[0].state.lp = 17.0;
                cells[0].lp_at_grace_end = 17.0;
            }
            observer.after_biochemistry(cells, step);
        }
        let mut damp = vec![0.0; cells.len()];
        damp[0] = 25.5;
        observer.before_terminal(cells, &damp, 3);
        damp[1] = 18.0;
        observer.after_terminal(&damp);
        let complete = &observer.ferroptotic_events[0];
        assert_eq!(complete.death_lp, 11.0);
        assert_eq!(complete.release_step, Some(2));
        assert_eq!(complete.release_lp, Some(17.0));
        assert_eq!(complete.release_damp, Some(25.5));
        assert_eq!(complete.horizon_lp, None);
        let censored = &observer.ferroptotic_events[1];
        assert_eq!(censored.scheduled_release_step, 3);
        assert_eq!(censored.release_lp, None);
        assert_eq!(censored.horizon_lp, Some(12.0));
        assert_eq!(censored.terminal_damp, Some(18.0));
        assert_eq!(observer.terminal.censored_deaths, 1);
        assert_eq!(observer.terminal.terminal_additions, 1);
        assert_eq!(observer.terminal.damp_before_terminal, 25.5);
        assert_eq!(observer.terminal.damp_after_terminal, 43.5);
    }

    #[test]
    fn empty_populations_remain_empty() {
        let grid = TumorGrid3D::generate(3, 3, 3, 20.0, 42);
        let cells = &grid.cells;
        let damp = vec![0.0; cells.len()];
        let mut observer = Measurements::new(
            cells.len(),
            5,
            1.0,
            IMMUNE_START_STEP,
            DAMP_KILL_THRESHOLD,
            None,
        );
        observer.after_biochemistry(cells, 0);
        observer.before_immune(cells, &damp, 0);
        observer.after_immune(cells, &damp, 0);
        observer.before_terminal(cells, &damp, 1);
        observer.after_terminal(&damp);
        assert!(observer.ferroptotic_events.is_empty());
        assert!(observer.eligible_cells.is_empty());
        assert!(observer.immune_kill_events.is_empty());
        assert_eq!(observer.terminal.terminal_damp, 0.0);
    }

    #[test]
    fn activation_uses_eligible_draws_and_inclusive_kd_boundaries() {
        let mut grid = TumorGrid3D::generate(3, 3, 3, 20.0, 42);
        for (idx, cell) in grid.cells.iter_mut().enumerate() {
            cell.is_tumor = idx < 6;
            cell.state.dead = false;
        }
        let mut damp = vec![0.0; grid.cells.len()];
        damp[..6].copy_from_slice(&[0.005, 0.01, 49.999, 50.0, 449.999, 450.0]);
        let mut observer = Measurements::new(grid.cells.len(), 5, 1.0, 1, 0.01, Some(50.0));
        for step in 0..3 {
            observer.after_biochemistry(&grid.cells, step);
            observer.before_immune(&grid.cells, &damp, step);
            observer.after_immune(&grid.cells, &damp, step);
        }
        observer.before_terminal(&grid.cells, &damp, 3);
        observer.after_terminal(&damp);
        let empty = observer.steps[0].activation.as_ref().unwrap();
        assert_eq!(empty.max_local_damp, None);
        assert_eq!(empty.activation_sum, 0.0);
        assert_eq!(empty.damp_ge_kd_opportunities, 0);
        assert_eq!(empty.damp_ge_9kd_opportunities, 0);
        let row = &observer.steps[1];
        let stats = row.activation.as_ref().unwrap();
        assert_eq!(row.eligible_cells, 5);
        assert_eq!(stats.max_local_damp, Some(450.0));
        assert_eq!(stats.damp_ge_kd_opportunities, 3);
        assert_eq!(stats.damp_ge_9kd_opportunities, 1);
        let expected: f64 = damp[1..6].iter().map(|d| d / (d + 50.0)).sum();
        assert_eq!(stats.activation_sum, expected);
        assert_eq!(observer.eligible_cells.len(), 5);
        for cell in &observer.eligible_cells {
            let stats = cell.activation.as_ref().unwrap();
            let local = damp[cell.cell_index];
            assert_eq!(cell.opportunities, 2);
            assert_eq!(stats.max_local_damp, Some(local));
            assert_eq!(stats.activation_sum, 2.0 * local / (local + 50.0));
            assert_eq!(
                stats.damp_ge_kd_opportunities,
                2 * usize::from(local >= 50.0)
            );
            assert_eq!(
                stats.damp_ge_9kd_opportunities,
                2 * usize::from(local >= 450.0)
            );
        }
    }

    #[test]
    fn disabled_activation_preserves_legacy_v1_json_bytes() {
        let grid = TumorGrid3D::generate(3, 3, 3, 20.0, 42);
        let damp = vec![0.0; grid.cells.len()];
        let mut observer = Measurements::new(grid.cells.len(), 5, 1.0, 60, 0.01, None);
        observer.after_biochemistry(&grid.cells, 0);
        observer.before_immune(&grid.cells, &damp, 0);
        observer.after_immune(&grid.cells, &damp, 0);
        observer.before_terminal(&grid.cells, &damp, 1);
        observer.after_terminal(&damp);
        assert_eq!(serde_json::to_string(&observer).unwrap(),
            "{\"ferroptotic_events\":[],\"immune_kill_events\":[],\"eligible_cells\":[],\"steps\":[{\"step\":0,\"ferroptotic_deaths\":0,\"completed_releases\":0,\"released_damp\":0.0,\"eligible_cells\":0,\"eligible_local_damp_sum\":0.0,\"immune_kills\":0}],\"terminal\":{\"censored_deaths\":0,\"terminal_additions\":0,\"terminal_damp\":-0.0,\"damp_before_terminal\":0.0,\"damp_after_terminal\":0.0}}");
        let cell = EligibleCell {
            cell_index: 3,
            first_step: 60,
            last_step: 61,
            opportunities: 2,
            local_damp_sum: 5.0,
            activation: None,
        };
        assert_eq!(serde_json::to_string(&cell).unwrap(),
            "{\"cell_index\":3,\"first_step\":60,\"last_step\":61,\"opportunities\":2,\"local_damp_sum\":5.0}");
    }

    #[test]
    fn repeated_opportunities_count_one_cell_and_immune_death_releases_nothing() {
        let mut grid = TumorGrid3D::generate(3, 3, 3, 20.0, 42);
        let cells = &mut grid.cells;
        for cell in cells.iter_mut() {
            cell.is_tumor = false;
        }
        cells[0].is_tumor = true;
        cells[0].state.dead = false;
        cells[0].state.death_step = None;
        let damp = vec![DAMP_KILL_THRESHOLD; cells.len()];
        let mut observer = Measurements::new(
            cells.len(),
            5,
            1.0,
            IMMUNE_START_STEP,
            DAMP_KILL_THRESHOLD,
            None,
        );
        let horizon = IMMUNE_START_STEP + 3;
        for step in 0..horizon {
            observer.after_biochemistry(cells, step);
            observer.before_immune(cells, &damp, step);
            if step == IMMUNE_START_STEP + 1 {
                // Exactly how the production immune loop records its death.
                cells[0].state.dead = true;
            }
            observer.after_immune(cells, &damp, step);
        }
        observer.before_terminal(cells, &damp, horizon);
        observer.after_terminal(&damp);
        assert_eq!(observer.eligible_cells.len(), 1);
        assert_eq!(observer.eligible_cells[0].opportunities, 2);
        assert_eq!(observer.eligible_cells[0].first_step, IMMUNE_START_STEP);
        assert_eq!(observer.eligible_cells[0].last_step, IMMUNE_START_STEP + 1);
        assert_eq!(observer.steps[IMMUNE_START_STEP as usize].eligible_cells, 1);
        assert_eq!(
            observer.steps[(IMMUNE_START_STEP + 1) as usize].immune_kills,
            1
        );
        assert_eq!(
            observer.steps[(IMMUNE_START_STEP + 2) as usize].eligible_cells,
            0
        );
        assert_eq!(observer.immune_kill_events.len(), 1);
        assert!(observer.ferroptotic_events.is_empty());
        assert_eq!(observer.terminal.terminal_additions, 0);
    }
}
