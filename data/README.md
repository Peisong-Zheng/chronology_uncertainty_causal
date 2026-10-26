# Local generated data

Windows is the authoritative computation platform from the 2026-10-26 migration onward. No Mac simulation inputs or chronology arrays are required to start a fresh Windows run. Generate them with `02_check_causal_pipeline.py`, then run `02_3_check_uncertainty_escalation.py` with the same ensemble size.

- `processed/<script_stem>/`: generated inputs, parameters, member results and summaries.
- `cache/`: native GPDC null distributions and their matching metadata. Regenerate these on a new platform.
- `logs/`: logs written by the formal launcher.
- `benchmarks/`: worker-count measurements, separate from scientific experiment outputs.
- `legacy_macos/`: preserved historical Mac outputs, including the stopped experiment05. These are not silently reused by the new Windows runs.

Generated data stay outside Git. Ignoring them does not delete local files. For a large completed run, save the full output folder separately if needed; its parameters, input arrays and member results belong together. Do not sync a folder while another machine is writing to it.

To start a different scientific configuration, preserve or rename the previous output folder before creating a fresh one. Resume is for the same configuration and environment, not for mixing experiments.
