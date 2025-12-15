use std::io::{self, Write};
use std::time::{Duration, Instant};

use indicatif::{ProgressBar, ProgressStyle};

pub struct StepSpinner {
    bar: Option<ProgressBar>,
    start: Instant,
    step: u8,
    total_steps: u8,
    step_start: Instant,
}

impl StepSpinner {
    pub fn new(total_steps: u8) -> Self {
        let now = Instant::now();
        Self {
            bar: None,
            start: now,
            step: 0,
            total_steps,
            step_start: now,
        }
    }

    pub fn step(&mut self, description: &str) {
        if let Some(bar) = self.bar.take() {
            bar.finish_and_clear();
        }

        self.step += 1;
        self.step_start = Instant::now();

        let bar = ProgressBar::new_spinner();
        bar.set_style(
            ProgressStyle::default_spinner()
                .template("  {spinner:.cyan} {msg}")
                .expect("invalid template")
                .tick_chars("⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"),
        );
        bar.enable_steady_tick(Duration::from_millis(80));
        bar.set_message(format!(
            "[{}/{}] {}...",
            self.step, self.total_steps, description
        ));

        self.bar = Some(bar);
    }

    pub fn complete_step(&mut self, description: &str, substeps: &[&str]) {
        if let Some(bar) = self.bar.take() {
            bar.finish_and_clear();
        }

        let elapsed = self.step_start.elapsed();
        let mut stderr = io::stderr().lock();

        let _ = writeln!(
            stderr,
            "  \x1b[32m✓\x1b[0m {:<44} {:>5.1}s",
            description,
            elapsed.as_secs_f64()
        );

        for substep in substeps {
            let _ = writeln!(stderr, "      \x1b[2m·\x1b[0m {}", substep);
        }
    }

    pub fn finish(mut self) {
        if let Some(bar) = self.bar.take() {
            bar.finish_and_clear();
        }

        print_footer(self.start.elapsed());
    }
}

fn print_footer(elapsed: Duration) {
    let mut stderr = io::stderr().lock();

    let _ = writeln!(stderr);
    let _ = writeln!(
        stderr,
        "  \x1b[2m╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╸\x1b[0m"
    );
    let _ = writeln!(stderr);
    let _ = writeln!(
        stderr,
        "  \x1b[32m✓\x1b[0m Forge complete {:>35}",
        format!("Total: {:.2}s", elapsed.as_secs_f64())
    );
    let _ = writeln!(stderr);
}

pub struct SilentProgress {}

impl SilentProgress {
    pub fn new() -> Self {
        Self {}
    }

    pub fn step(&mut self, _description: &str) {}

    pub fn complete_step(&mut self, _description: &str, _substeps: &[&str]) {}
}

impl Default for SilentProgress {
    fn default() -> Self {
        Self::new()
    }
}

pub enum Progress {
    Interactive(StepSpinner),
    Silent(SilentProgress),
}

impl Progress {
    pub fn new(interactive: bool, total_steps: u8) -> Self {
        if interactive {
            Self::Interactive(StepSpinner::new(total_steps))
        } else {
            Self::Silent(SilentProgress::new())
        }
    }

    pub fn step(&mut self, description: &str) {
        match self {
            Self::Interactive(s) => s.step(description),
            Self::Silent(s) => s.step(description),
        }
    }

    pub fn complete_step(&mut self, description: &str, substeps: &[&str]) {
        match self {
            Self::Interactive(s) => s.complete_step(description, substeps),
            Self::Silent(s) => s.complete_step(description, substeps),
        }
    }

    pub fn finish(self) {
        match self {
            Self::Interactive(s) => s.finish(),
            Self::Silent(_) => {}
        }
    }
}
