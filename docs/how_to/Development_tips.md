# Development tips

This page describes any tips we have for developers.

## Cppcheck

Our CI uses [cppcheck](https://cppcheck.sourceforge.io/) to verify that best practices for C++ are followed. In addition we have some addons that have been written to catch common coding errors. To run these tests for Gyselalib++ or for a project which has Gyselalib++ as a submodule, simply use the script `${GYSELALIB_ROOT}/bin/run_cppcheck`.
