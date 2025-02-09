import logging
from pathlib import Path
from unittest.mock import MagicMock, patch

from legenddataflow.pars_loading import ParsCatalog

log = logging.getLogger(__name__)


def test_match_pars_files():
    filelist1 = ["file1_step1_typeA", "file2_step2_typeB"]
    filelist2 = ["file3_step1_typeA", "file4_step3_typeC"]

    with patch(
        "legenddataflow.pars_loading.ProcessingFileKey.get_filekey_from_pattern"
    ) as mock_get_filekey:
        mock_get_filekey.side_effect = [
            MagicMock(processing_step="step1", datatype="typeA"),
            MagicMock(processing_step="step2", datatype="typeB"),
            MagicMock(processing_step="step1", datatype="typeA"),
            MagicMock(processing_step="step3", datatype="typeC"),
        ]

        result1, result2 = ParsCatalog.match_pars_files(filelist1, filelist2)

        assert result1 == ["file1_step1_typeA", "file3_step1_typeA"]
        assert result2 == ["file4_step3_typeC"]


def test_get_par_file():
    setup = "test_setup"
    timestamp = "20230101T000000Z"
    tier = "test_tier"

    catalog = ParsCatalog.get(
        [
            {"valid_from": "20230101T000000Z", "apply": ["file1.yaml", "file2.yaml"]},
        ]
    )
    log.debug(catalog)

    with patch(
        "legenddataflow.pars_loading.get_pars_path"
    ) as mock_get_pars_path, patch(
        "legenddataflow.pars_loading.par_overwrite_path"
    ) as mock_par_overwrite_path, patch(
        "legenddataflow.pars_loading.ParsCatalog.get_files"
    ) as mock_get_files, patch(
        "legenddataflow.pars_loading.ParsCatalog.match_pars_files"
    ) as mock_match_pars_files:
        mock_get_pars_path.return_value = "/pars/path"
        mock_par_overwrite_path.return_value = "/overwrite/path"
        mock_get_files.side_effect = [
            ["file1.yaml", "file2.yaml"],  # pars_files
            ["file3.yaml"],  # pars_files_overwrite
        ]
        mock_match_pars_files.return_value = (
            ["file1.yaml", "file2.yaml"],
            ["file3.yaml"],
        )

        result = ParsCatalog.get_par_file(catalog, setup, timestamp, tier)

        expected_result = [
            Path("/pars/path/file1.yaml"),
            Path("/pars/path/file2.yaml"),
            Path("/overwrite/path/test_tier/file3.yaml"),
        ]

        assert result == expected_result
