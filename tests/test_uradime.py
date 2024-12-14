import pytest
from Bio.Seq import Seq
from URAdime.URAdime import (
    load_primers,
    is_match,
    check_terminal_match,
    get_base_primer_name,
    format_summary_table
)
import pandas as pd
import os

# Test data
@pytest.fixture
def sample_primers_df():
    data = {
        'Name': ['Primer1', 'Primer2'],
        'Forward': ['ACGTACGT', 'GCTAGCTA'],
        'Reverse': ['TGCATGCA', 'CGATCGAT'],
        'Size': [100, 200]
    }
    return pd.DataFrame(data)

def test_is_match():
    assert is_match("ACGTACGT", "ACGTACGT", 0) == True
    assert is_match("ACGTACGT", "ACGTACGA", 1) == True
    assert is_match("ACGTACGT", "ACGTACGA", 0) == False
    assert is_match("", "ACGTACGT", 0) == False
    assert is_match("ACGTACGT", "", 0) == False

def test_check_terminal_match():
    seq = "ACGTACGTNNNN"
    primer = "ACGTACGT"
    found, length = check_terminal_match(seq, primer)
    assert found == True
    assert length == 8

def test_get_base_primer_name():
    assert get_base_primer_name("Primer1_Forward") == "Primer1"
    assert get_base_primer_name("Primer2_Reverse_Terminal_15bp") == "Primer2"
    assert get_base_primer_name("None") == None

def test_format_summary_table():
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