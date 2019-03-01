from EnteroTyper.bin.accessories import run_subprocess, dependency_check, check_all_dependencies


def test_run_subprocess():
    assert run_subprocess('echo test', get_stdout=True) == 'test'
    assert run_subprocess('echo test', get_stdout=False) is None


def test_dependency_check():
    assert dependency_check('echo') is True
    assert dependency_check('ls') is True
    assert dependency_check('DefinitelyNotARealProgram__123') is False


def test_check_all_dependencies():
    assert check_all_dependencies(['echo', 'ls']) is True
    assert check_all_dependencies(['Definitely__Fake1', 'Definitely__Fake2']) is False
