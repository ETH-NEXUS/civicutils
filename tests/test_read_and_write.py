import os
import pytest
from civicutils.read_and_write import read_in_snvs, read_in_cnvs, read_in_expr

# Define paths to example data files
snv_file = os.path.join("civicutils","data", "example_snv.txt")
cnv_file = os.path.join("civicutils","data", "example_cnv.txt")
expr_file = os.path.join("civicutils","data", "example_expr.txt")


def test_read_in_snvs():
    raw_data, snv_data, extra_header = read_in_snvs(snv_file)

    # Basic assertions to ensure data is read correctly
    assert isinstance(raw_data, dict)
    assert isinstance(snv_data, dict)
    assert isinstance(extra_header, list)

    # Check if expected keys exist
    assert 'ATR' in snv_data  
    assert 'c.5488C>T|p.Gln1830Ter|||1' in snv_data['ATR']

def test_read_in_cnvs():
    raw_data, cnv_data, extra_header = read_in_cnvs(cnv_file)

    # Basic assertions to ensure data is read correctly
    assert isinstance(raw_data, dict)
    assert isinstance(cnv_data, dict)
    assert isinstance(extra_header, list)

    # Check if expected keys exist
    assert 'BRCA1' in cnv_data  # Replace with actual gene names in your data
    assert 'AMPLIFICATION|0' in cnv_data['BRCA1']

def test_read_in_expr():
    raw_data, expr_data, extra_header = read_in_expr(expr_file)

    # Basic assertions to ensure data is read correctly
    assert isinstance(raw_data, dict)
    assert isinstance(expr_data, dict)
    assert isinstance(extra_header, list)

    # Check if expected keys exist
    assert 'TP53' in expr_data  # Replace with actual gene names in your data
    assert '-2.32156489488682|4' in expr_data['TP53']

if __name__ == "__main__":
    pytest.main()