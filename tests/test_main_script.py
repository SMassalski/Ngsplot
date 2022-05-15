import os
import pytest
from make_plots import Config, get_cli_args


# TODO: Keeps defaults test

data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data'))


@pytest.fixture
def config():
    return Config()


def test_get_cli_args_bool(config):
    args = get_cli_args(['--save_data', '--sort', '--only_first'])
    config.update(vars(args))
    assert config.save_data is True
    assert config.sort is True
    assert config.only_first is True


# Config.convert_type()
def test_config_convert_true(config):
    assert config.convert_type('true') is True
    assert config.convert_type('True') is True
    assert config.convert_type('TRUE') is True
    
    
def test_config_convert_false(config):
    assert config.convert_type('false') is False
    assert config.convert_type('False') is False
    assert config.convert_type('FALSE') is False


def test_config_convert_int(config):
    assert config.convert_type('1') == 1


def test_config_convert_float(config):
    assert config.convert_type('1.25') == 1.25
    
    
def test_config_convert_none(config):
    assert config.convert_type('none') is None
    assert config.convert_type('None') is None
    assert config.convert_type('NONE') is None


def test_config_convert_list(config):
    result = config.convert_type('item1,item2,item3')
    assert result == ['item1', 'item2', 'item3']
    
    
def test_config_convert_str(config):
    assert config.convert_type('string') == 'string'


# Config.update()
def test_config_update_none(config):
    config.flank = 3000
    config.update(dict(flank=None))
    assert config.flank is not None
    
    
def test_config_update_true(config):
    config.update(dict(sort=True))
    assert config.sort is True
    

def test_config_update_false(config):
    config.update(dict(sort=False))
    assert config.sort is False
    config.sort = True
    config.update(dict(sort=False))
    assert config.sort is True


def test_config_update_str(config):
    config.update(dict(out_format='pdf'))
    assert config.out_format == 'pdf'


def test_config_update_smooth_none(config):
    config.update(dict(smooth='none'))
    assert config.smooth is None


def test_config_update_smooth_true(config):
    config.update(dict(smooth='true'))
    assert config.smooth is True


def test_config_update_smooth_float(config):
    config.update(dict(smooth='0.12'))
    assert config.smooth == 0.12
    

@pytest.fixture
def from_file(config, request, caplog):
    marker = request.node.get_closest_marker('fp')
    if marker is None:
        fp = 'correct_test_config.txt'
    else:
        fp = marker.args[0]
    config.from_file(os.path.join(data_dir, fp))
    return {'logs': caplog.text}


def test_config_from_file_correct_str(config, from_file):
    assert config.signal_file == 'fake_file.bedGraph'
    assert config.roi == 'end'
    
    
def test_config_from_file_correct_bool(config, from_file):
    assert config.sort is True
    assert config.colorbar is False


def test_config_from_file_correct_int(config, from_file):
    assert isinstance(config.flank, int) and config.flank == 1


def test_config_from_file_correct_float(config, from_file):
    assert isinstance(config.min_score, float) and config.min_score == 10.13
    

def test_config_from_file_correct_list(config, from_file):
    assert config.omit_chr == ['chr1', 'chr2']


@pytest.mark.fp('bad_config.txt')
def test_config_from_file_missing_line(from_file):
    assert 'Wrong number of columns in line 3' in from_file['logs']


@pytest.mark.fp('bad_config.txt')
def test_config_from_file_no_value(from_file):
    assert 'Wrong number of columns in line 2' in from_file['logs']
    
    
@pytest.mark.fp('bad_config.txt')
def test_config_from_file_to_many_values(from_file):
    assert 'Wrong number of columns in line 6' in from_file['logs']


@pytest.mark.fp('bad_config.txt')
def test_config_from_file_unknown_argument(from_file):
    assert 'Unrecognized argument "fake_argument" in line 5' \
           in from_file['logs']
