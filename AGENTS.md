# AGENTS.md

## Critical Rules

1. Add unit tests for any new functionality you implement. Ensure that the tests cover various edge cases and scenarios to validate the correctness of your code.
2. The project uses `pixi` for package management. Use the tool for building (`pixi run build`), testing (`pixi run test`) and formatting (`pixi run format`).
3. Always format the code.
4. Ensure all tests pass before finishing your work. Run `pixi run test` to verify that your changes do not break existing functionality.
5. Never commit any code.
6. Delete any temporary files generated while implementing new functionality.

## Preferences

1. Use `nalgebra` for mathematical operations and data structures. Avoid using other libraries for these purposes.
