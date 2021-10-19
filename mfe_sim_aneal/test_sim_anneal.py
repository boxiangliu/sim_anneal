import pytest
from sim_anneal import *

CODING_WHEEL_FN = "coding_wheel.txt"
PROTEIN_FN = "membrane.protein"
CODON_TABLE = read_coding_wheel(CODING_WHEEL_FN)

def test_read_coding_wheel():
	codon_table = read_coding_wheel(CODING_WHEEL_FN)
	assert codon_table["Phe"] == ["UUU", "UUC"]
	assert codon_table["STOP"] == ["UAA", "UAG", "UGA"]

def test_protein():
	protein = read_protein(PROTEIN_FN)
	assert protein[0] == "Met"
	assert protein[-1] == "STOP"
	assert protein[1] == "Ala"

def test_initialize():
	protein = ["Met", "Trp"]
	rna = RNA(protein, CODON_TABLE)
	assert rna.rna == ["AUG", "UGG"]

def test_mfe():
	rna = RNA("", CODON_TABLE)
	rna.rna = "CCCCCCCCCGGGGGGGGG"
	rna.get_mfe("test.txt")
	assert rna.mfe == -15.40

def test_mutate():
	rna = RNA(["Met", "Phe"], CODON_TABLE)
	rna.rna = ['AUG', 'UUU']
	rna.mutate()
	assert rna.rna == ['AUG', 'UUC']