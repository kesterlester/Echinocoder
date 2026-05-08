def pytest_report_header(config):
    return (
        "NOTE: MinimalConfusableSets/ is excluded from top-level pytest collection "
        "because its bistate/tristate config flags conflict when both test modules are "
        "loaded in the same process. Run its tests via MinimalConfusableSets/run_tests.sh "
        "instead. To re-enable top-level recursion, remove 'MinimalConfusableSets' from "
        "norecursedirs in pytest.ini."
    )
