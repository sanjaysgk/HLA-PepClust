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
    (res_dir / "test_3g.ds").write_text(
        """HEADER
pos A C D E
1 0.1 0.2 0.3 0.4
2 0.5 0.6 0.7 0.8
""", encoding='utf-8'
    )
    return str(temp_dir)


@pytest.fixture
def hla_reference_files(temp_dir):
    # Create fake HLA reference files
    ref_dir = temp_dir / "human_ref"
    ref_dir.mkdir()
    (ref_dir / "HLA_A0201.txt").write_text("dummy", encoding='utf-8')
    (ref_dir / "HLA_B0702.txt").write_text("dummy", encoding='utf-8')
    return str(ref_dir)

@pytest.fixture
def gibbs_files(temp_dir):
    # Create fake Gibbs output structure
    res_dir = temp_dir / "res"
    res_dir.mkdir(parents=True, exist_ok=True)
    
    # Create a fake gibbs.3g.ds.out file with sequence data and metadata
    gibbs_file_content = """# /tools/src/gibbscluster-2.0/bin/gibbs_cluster_1.2.2_DB.Linux_x86_64 -mhc1 -trash -l 9 -g 1 -ts 1.5 -blf /tools/src/gibbscluster-2.0/data/blosum%i.freq_rownorm -blm 62 -base 2 -bg 1 -i 10 -s 2 -nt 20 -fr 20 -fs 100 -wlc 200 -swt 0 -lambda 0.8 -dlen 4 -ilen 1 -findels 10 -sigma 5 /var/www/services/services/GibbsCluster-2.0/tmp/D90_HLA_3844874/data/pepdata.txt
    # Mon Feb 10 01:33:36 2025
    # User: services
    # PWD : /var/www/webface/tmp/server/GibbsCluster/67A9495E0012FB421004CC0F
    # Host: Linux engine 6.8.0-49-generic x86_64
    # Number of temperature steps 20 dt 0.078947
    # Background frequencies
    # 0 A  0.07400
    # 1 R  0.05200
    # 2 N  0.04500
    # 3 D  0.05400
    # ...
    pos  A     C     D     E
    1    0.1   0.2   0.3   0.4
    2    0.5   0.6   0.7   0.8
    3    0.9   1.0   1.1   1.2
    """
    
    # Create the gibbs.3g.ds.out file in the `res` directory
    file_path = res_dir / "gibbs.3g.ds.out"
    file_path.write_text(gibbs_file_content, encoding="utf-8")
    
    return str(temp_dir)  # Return the base directory for later use

def test_parse_gibbs_output_valid(cluster_search, gibbs_files):
    # Path to the fake Gibbs file for 3 clusters
    test_file = os.path.join(gibbs_files, "res", "gibbs.3g.ds.out")
    
    # Call the function to parse the Gibbs output
    df = cluster_search.parse_gibbs_output(test_file, 3)
    
    # Assertions to validate the output
    assert isinstance(df, pd.DataFrame)  # Ensure it is a DataFrame
    assert len(df) == 3  # Check that there are 3 rows (matching the 3 positions in the test file)
    assert list(df.columns) == ["pos", "A", "C", "D", "E"]  # Ensure the columns match

def test_check_HLA_DB_valid(cluster_search, hla_reference_files):
    assert cluster_search.check_HLA_DB(["HLA_A0201"], hla_reference_files) is True
    assert "HLA_A0201" in cluster_search.valid_HLA


def test_check_HLA_DB_invalid(cluster_search, hla_reference_files):
    assert cluster_search.check_HLA_DB(["HLA_INVALID"], hla_reference_files) is False


@patch.object(ClusterSearch, "_compute_and_log_correlation")
def test_compute_correlations_all(
    mock_compute, cluster_search, temp_dir, hla_reference_files
):
    # Setup test matrices
    gibbs_mat_dir = temp_dir / "matrices"
    gibbs_mat_dir.mkdir()
    (gibbs_mat_dir / "test1.mat").write_text("dummy", encoding='utf-8')
    (gibbs_mat_dir / "test2.mat").write_text("dummy", encoding='utf-8')

    cluster_search.compute_correlations(
        gibbs_folder=str(temp_dir),
        human_reference_folder=hla_reference_files,
        n_clusters="all",
        output_path=str(temp_dir),
        hla_list=None,
    )

    assert mock_compute.call_count >= 4  # 2 test mats * 2 ref mats


def test_format_input_gibbs(cluster_search, temp_dir):
    # Create the fake Gibbs matrix file (in the format you provided)
    test_file = temp_dir / "test.mat"
    test_file.write_text(
        """A R N D C Q E G H I L K M F P S T W Y V
1 L -0.003 -0.754 -0.565 0.534 -5.654 0.176 -2.475 -2.73 0.813 0.382 1.1 -2.714 0.255 0.66 -5.08 1.465 0.8 -3.204 1.528 0.937
2 A 2.192 -7.043 -3.805 -5.099 -9.259 3.42 -1.504 -0.112 -4.502 -1.436 -0.469 -8.574 1.038 -2.592 3.147 1.618 1.118 -4.424 -4.034 -1.697
3 L 0.293 -1.474 0.84 2.322 -9.317 -2.533 -0.221 -5.167 -2.0 0.613 1.068 -5.163 -0.312 2.681 -1.49 -2.173 -2.634 -0.02 3.177 0.283
4 P -1.285 -3.291 -1.228 0.33 -5.402 0.772 1.858 1.734 -0.895 -3.064 -2.082 -1.132 -2.816 -3.546 3.797 0.957 1.008 -1.586 -2.873 -3.702
5 L -0.095 -1.62 -0.444 -1.934 -7.538 1.413 -0.81 0.181 0.415 0.235 1.462 -4.064 -1.029 0.292 0.141 0.131 0.291 -2.087 -0.949 1.363
6 P 1.446 -3.604 -4.032 -4.336 -6.374 -0.749 -3.261 0.094 -2.445 1.838 -0.494 -5.305 0.693 1.462 3.553 -1.9 -0.013 -2.736 -1.629 1.835
7 L -1.269 -0.214 0.272 -0.638 -6.146 2.619 0.968 -6.535 -1.132 -1.61 2.328 -3.365 0.966 -2.825 -1.871 1.532 0.604 -2.84 -2.259 -0.708
8 T -0.016 -0.948 -1.943 -4.893 -4.627 0.82 1.449 -3.573 -0.252 -1.883 0.332 -0.971 -1.268 -0.04 1.944 0.22 2.179 -4.429 -0.213 0.97
9 V -2.663 -7.59 -8.374 -9.28 -10.726 -7.079 -10.786 -10.518 -6.559 4.128 1.621 -8.698 -1.603 1.725 -5.828 -7.278 -2.865 -7.609 3.969 4.372
""",
        encoding="utf-8",
    )
    
    # Test the formatting function
    df = cluster_search.format_input_gibbs(str(test_file))
    
    # Assert that the shape of the returned DataFrame is correct
    assert df.shape == (9, 20)  # There are 9 rows and 20 amino acids
    # Assert that the columns correspond to the amino acids in the first row
    expected_columns = [
        "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"
    ]
    assert list(df.columns) == expected_columns

    # Verify that the numeric conversion worked as expected
    assert df.iloc[0, 0] == -0.003  # Check the first numeric value
    assert df.iloc[1, 1] == -7.043  # Check another numeric value



@patch("PIL.Image.open")
@patch("os.listdir")
def test_create_image_grid(mock_listdir, mock_open, cluster_search, temp_dir):
    # Mock image data
    mock_img = MagicMock(spec=Image.Image)
    mock_img.width = 100
    mock_img.height = 100
    mock_open.return_value = mock_img

    # Mock directory listings
    mock_listdir.side_effect = [
        ["cluster1.mat.png"],  # gibbs images
        ["HLA_A0201.png"],  # reference images
    ]

    # Create dummy correlation data
    cluster_search.correlation_dict = {
        ("cluster1.mat", "HLA_A0201"): 0.95,
        ("cluster2.mat", "HLA_B0702"): 0.85,
    }

    cluster_search.valid_HLA = ["A0201"]

    cluster_search.create_image_grid(
        correlation_dict=cluster_search.correlation_dict,
        image_folder=str(temp_dir),
        DB_image_folder=str(temp_dir),
        output_path=str(temp_dir),
        HLA_list=["A0201"],
    )

    # Verify output file was created
    assert os.path.exists(
        os.path.join(str(temp_dir), "compare_allotypes.png")  # Corrected typo here
    )


def test_generate_unique_random_ids(cluster_search):
    ids = cluster_search.generate_unique_random_ids(10)
    assert len(ids) == 10
    assert len(set(ids)) == 10
    assert all(100000 <= id <= 999999 for id in ids)
