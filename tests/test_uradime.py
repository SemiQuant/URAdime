# tests/test_uradime.py
import pytest
import pandas as pd
import os
from URAdime.URAdime import (
    is_match,
    check_terminal_match,
    get_base_primer_name,
    format_summary_table
)

# Test fixtures
@pytest.fixture
def sample_primers_df():
    """Create a sample primers DataFrame for testing"""
    data = {
        'Name': ['Primer1', 'Primer2'],
        'Forward': ['ACGTACGT', 'GCTAGCTA'],
        'Reverse': ['TGCATGCA', 'CGATCGAT'],
        'Size': [100, 200]
    }
    return pd.DataFrame(data)

# Test functions
def test_is_match():
    """Test sequence matching function"""
    assert is_match("ACGTACGT", "ACGTACGT", 0) == True
    assert is_match("ACGTACGT", "ACGTACGA", 1) == True
    assert is_match("ACGTACGT", "ACGTACGA", 0) == False
    assert is_match("", "ACGTACGT", 0) == False
    assert is_match("ACGTACGT", "", 0) == False
    assert is_match("ACGTACGT", "NNNNNNNN", 0) == True

def test_check_terminal_match():
    """Test terminal sequence matching"""
    sequence = "ACGTACGTTTTT"
    primer = "ACGTACGT"
    found, length = check_terminal_match(sequence, primer, terminus_length=8)
    assert found == True
    assert length == 8

def test_get_base_primer_name():
    """Test primer name extraction"""
    assert get_base_primer_name("Primer1_Forward") == "Primer1"
    assert get_base_primer_name("Primer2_Reverse_Terminal_15bp") == "Primer2"
    assert get_base_primer_name("None") == None

def test_format_summary_table():
    """Test summary table formatting"""
    data = {
        'Category': ['Category1', 'Category2'],
        'Count': [10, 20],
        'Percentage': [25.5, 74.5]
    }
    df = pd.DataFrame(data)
    table = format_summary_table(df)
    assert isinstance(table, str)
    assert 'Category1' in table
    assert 'Category2' in table