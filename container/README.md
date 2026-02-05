# PANOSETI Docker Container

Minimal usage via the Makefile.

## Prereqs

- Docker installed and running
- CORSIKA credentials

## Build

```bash
export CORSIKA_USER=your_username
export CORSIKA_PASSWORD=your_password
make build
```

## Test

```bash
make test
```

## Shell

```bash
make run
```

## Example simulation (batch)

```bash
make example-sim
```

## Example visualization (requires X11)

```bash
make example-vis
```

## One-step example (simulation + visualization, requires X11)

```bash
make example
```

## Help

```bash
make help
```
