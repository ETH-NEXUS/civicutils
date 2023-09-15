from civicpy import civic
import pytest 

civic.load_cache("/cluster/customapps/biomed/nexus/sharedutils/db_cache/civicDB/cache.pkl")

@pytest.fixture
def all_genes():
    all_results = civic.get_all_genes()
    return all_results

def test_loading_cache():
    assert civic.load_cache("/cluster/customapps/biomed/nexus/sharedutils/db_cache/civicDB/cache.pkl") == True
    
def test_genes_attributs(all_genes):
    assert hasattr(all_genes[0], "id") == True
    assert hasattr(all_genes[0], "entrez_id") == True
    assert hasattr(all_genes[0], "aliases") == True
    assert hasattr(all_genes[0], "variants") == True

def test_variant_attributs(all_genes):
    assert hasattr(all_genes[0].variants[0], "molecular_profiles") == True
    assert hasattr(all_genes[0].variants[0], "id") == True
    assert hasattr(all_genes[0].variants[0], "name") == True
    assert hasattr(all_genes[0].variants[0], "hgvs_expressions") == True
    assert hasattr(all_genes[0].variants[0], "variant_types") == True

def test_molecular_profil_attributs(all_genes):
    assert hasattr(all_genes[0].variants[0].molecular_profiles[0], "evidence_items") == True
    assert hasattr(all_genes[0].variants[0].molecular_profiles[0], "id") == True
    assert hasattr(all_genes[0].variants[0].molecular_profiles[0], "name") == True
    assert hasattr(all_genes[0].variants[0].molecular_profiles[0], "molecular_profile_score") == True

def test_evidence_items_attributs(all_genes):
    assert hasattr(all_genes[0].variants[0].molecular_profiles[0].evidence_items[0], "variant_origin") == True
    assert hasattr(all_genes[0].variants[0].molecular_profiles[0].evidence_items[0], "evidence_type") == True
    assert hasattr(all_genes[0].variants[0].molecular_profiles[0].evidence_items[0], "status") == True
    assert hasattr(all_genes[0].variants[0].molecular_profiles[0].evidence_items[0].source, "source_type") == True
    assert hasattr(all_genes[0].variants[0].molecular_profiles[0].evidence_items[0].source, "citation_id") == True
    assert hasattr(all_genes[0].variants[0].molecular_profiles[0].evidence_items[0], "therapies") == True
    
    
    
    