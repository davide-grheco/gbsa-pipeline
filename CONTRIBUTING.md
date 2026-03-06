# Contributing

Contributions are welcome, and they are greatly appreciated! Every little bit helps, and credit will always be given.

## Environment setup

Nothing easier!

Fork and clone the repository, then:

```bash
cd gbsa-pipeline
pixi install -e dev
pixi run -e dev pre-commit-install
```

You now have the dependencies installed in a pixi environment. Common commands:

- Run the app: `pixi run -e dev gbsa-pipeline [ARGS...]`
- Format: `pixi run -e dev format`
- Checks: `pixi run -e dev check`
- Tests: `pixi run -e dev test`
- Full CI set: `pixi run -e dev ci`
- Pre-commit (all files): `pixi run -e dev pre-commit-run`

## Tasks

Commands and tasks are exposed through pixi (see `pyproject.toml` under `tool.pixi.tasks`).

## Development

As usual:

1. create a new branch: `git switch -c feature-or-bugfix-name`
1. edit the code and/or the documentation

**Before committing:**

1. run `pixi run -e dev format` to auto-format the code
1. run `pixi run -e dev check` to check everything (fix any warning)
1. run `pixi run -e dev test` to run the tests (fix any issue)
1. if you updated the documentation or the project dependencies:
    1. run `pixi run -e dev docs`
    1. go to http://localhost:8000 and check that everything looks good
