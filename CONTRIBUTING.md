# Contributing to MPC Spacecraft Control

Thank you for your interest in contributing to this project! We welcome contributions from the community.

## How to Contribute

### Reporting Bugs

If you find a bug, please create an issue with:
- A clear, descriptive title
- Steps to reproduce the issue
- Expected behavior
- Actual behavior
- MATLAB version and relevant toolbox versions
- Any relevant log output or error messages

### Suggesting Enhancements

Enhancement suggestions are tracked as GitHub issues. When creating an enhancement suggestion, please include:
- A clear, descriptive title
- A detailed description of the proposed enhancement
- Explain why this enhancement would be useful
- List any alternative solutions you've considered

### Pull Requests

1. Fork the repository
2. Create a new branch (`git checkout -b feature/amazing-feature`)
3. Make your changes
4. Test your changes thoroughly with MATLAB
5. Commit your changes (`git commit -m 'Add some amazing feature'`)
6. Push to the branch (`git push origin feature/amazing-feature`)
7. Open a Pull Request

#### Pull Request Guidelines

- Follow the existing code style and conventions
- Add comments to your code where necessary
- Update documentation if you're changing functionality
- Test your changes with different parameter sets
- Keep pull requests focused - one feature/fix per PR

## Code Style

- Use meaningful variable names
- Follow MATLAB naming conventions:
  - Functions: `lowercase_with_underscores`
  - Variables: `camelCase` or `snake_case` (be consistent)
  - Constants: `UPPERCASE_WITH_UNDERSCORES`
- Add function headers with clear descriptions of inputs and outputs
- Comment complex algorithms or non-obvious code

## Testing

Before submitting a pull request:
1. Run `main_run_sim.m` with both baseline and MPC controllers
2. Verify the Simulink model builds and runs: `build_simulink_model()`
3. Test with different parameter configurations
4. Check that plots are generated correctly

## Questions?

Feel free to open an issue with the "question" label if you have any questions about contributing.

## Code of Conduct

This project follows a Code of Conduct that all contributors are expected to adhere to. Please read [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md) before contributing.
