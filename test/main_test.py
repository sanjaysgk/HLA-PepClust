import pytest
import os
import pandas as pd
import numpy as np
from unittest.mock import patch, MagicMock
from PIL import Image
from cli.cluster_search import ClusterSearch

@pytest.fixture
def cluster_search():
    return ClusterSearch()

@pytest.fixture
def temp_dir(tmp_path):
    d = tmp_path / "test_data"
    d.mkdir()
    return d

@pytest.fixture
def gibbs_files(temp_dir):
    # Create fake Gibbs output structure
    res_dir = temp_dir / "res"
    res_dir.mkdir()
    (res_dir / "test_3g.ds").write_text("""HEADER
pos A C D E
1 0.1 0.2 0.3 0.4
2 0.5 0.6 0.7 0.8
""")
    return str(temp_dir)

@pytest.fixture
def hla_reference_files(temp_dir):
    # Create fake HLA reference files
    ref_dir = temp_dir / "human_ref"
    ref_dir.mkdir()
    (ref_dir / "HLA_A0201.txt").write_text("dummy")
    (ref_dir / "HLA_B0702.txt").write_text("dummy")
    return str(ref_dir)

def test_parse_gibbs_output_valid(cluster_search, gibbs_files):
    df = cluster_search.parse_gibbs_output(gibbs_files, 3)
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 2
    assert list(df.columns) == ['pos', 'A', 'C', 'D', 'E']

def test_check_HLA_DB_valid(cluster_search, hla_reference_files):
    assert cluster_search.check_HLA_DB(["HLA_A0201"], hla_reference_files) is True
    assert "HLA_A0201" in cluster_search.valid_HLA

def test_check_HLA_DB_invalid(cluster_search, hla_reference_files):
    assert cluster_search.check_HLA_DB(["HLA_INVALID"], hla_reference_files) is False

@patch.object(ClusterSearch, '_compute_and_log_correlation')
def test_compute_correlations_all(mock_compute, cluster_search, temp_dir, hla_reference_files):
    # Setup test matrices
    gibbs_mat_dir = temp_dir / "matrices"
    gibbs_mat_dir.mkdir()
    (gibbs_mat_dir / "test1.mat").write_text("dummy")
    (gibbs_mat_dir / "test2.mat").write_text("dummy")
    
    cluster_search.compute_correlations(
        gibbs_folder=str(temp_dir),
        human_reference_folder=hla_reference_files,
        n_clusters="all",
        output_path=str(temp_dir),
        hla_list=None
    )
    
    assert mock_compute.call_count >= 4  # 2 test mats * 2 ref mats

def test_format_input_gibbs(cluster_search, temp_dir):
    test_file = temp_dir / "test.mat"
    test_file.write_text("A C D E\n0.1 0.2 0.3 0.4\n0.5 0.6 0.7 0.8")
    df = cluster_search.format_input_gibbs(str(test_file))
    assert df.shape == (2, 4)
    assert list(df.columns) == ['A', 'C', 'D', 'E']

@patch('PIL.Image.open')
@patch('os.listdir')
def test_create_image_grid(mock_listdir, mock_open, cluster_search, temp_dir):
    # Mock image data
    mock_img = MagicMock(spec=Image.Image)
    mock_img.width = 100
    mock_img.height = 100
    mock_open.return_value = mock_img
    
    # Mock directory listings
    mock_listdir.side_effect = [
        ['cluster1.mat.png'],  # gibbs images
        ['HLA_A0201.png']      # reference images
    ]
    
    # Create dummy correlation data
    cluster_search.correlation_dict = {
        ('cluster1.mat', 'HLA_A0201'): 0.95,
        ('cluster2.mat', 'HLA_B0702'): 0.85
    }
    
    cluster_search.valid_HLA = ['A0201']
    
    cluster_search.create_image_grid(
        correlation_dict=cluster_search.correlation_dict,
        image_folder=str(temp_dir),
        DB_image_folder=str(temp_dir),
        output_path=str(temp_dir),
        HLA_list=['A0201']
    )
    
    # Verify output file was created
    assert os.path.exists(os.path.join(cluster_search.full_path, "campare_allotypes.png"))

def test_generate_unique_random_ids(cluster_search):
    ids = cluster_search.generate_unique_random_ids(10)
    assert len(ids) == 10
    assert len(set(ids)) == 10
    assert all(100000 <= id <= 999999 for id in ids)