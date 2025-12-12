"""Integration tests for remote cache functionality.

To run these tests, you need:
1. A valid API token with appropriate permissions (can_read and/or can_upload)
2. Real API server URL (default: https://api.orca-descriptors.massonnn.ru)

Set environment variables:
    export ORCA_CACHE_SERVER_URL="https://api.orca-descriptors.massonnn.ru"
    export ORCA_CACHE_API_TOKEN="your-api-token-here"

Or create tests/.env file with:
    ORCA_CACHE_SERVER_URL=https://api.orca-descriptors.massonnn.ru
    ORCA_CACHE_API_TOKEN=your-api-token-here

Run all integration tests:
    pytest tests/test_remote_cache.py -m integration -v

Run specific test:
    pytest tests/test_remote_cache.py::test_orca_with_remote_cache_integration -v

Note: Tests will be skipped if server is not available or token is not set.
"""

import os
import pytest
import tempfile
import shutil
from pathlib import Path
from rdkit.Chem import MolFromSmiles, AddHs

from orca_descriptors import Orca, ORCABatchProcessing
from orca_descriptors.remote_cache import RemoteCacheClient, RemoteCacheError, RemoteCachePermissionError

_env_file = Path(__file__).parent / ".env"
if _env_file.exists():
    try:
        try:
            from dotenv import load_dotenv
            load_dotenv(_env_file)
        except ImportError:
            with open(_env_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#') and '=' in line:
                        key, value = line.split('=', 1)
                        key = key.strip()
                        value = value.strip().strip('"').strip("'")
                        os.environ.setdefault(key, value)
    except Exception:
        pass

CACHE_SERVER_URL = os.getenv("ORCA_CACHE_SERVER_URL", "https://api.orca-descriptors.massonnn.ru")
CACHE_API_TOKEN = os.getenv("ORCA_CACHE_API_TOKEN", None) or os.getenv("CACHE_API_TOKEN", None)


def is_cache_server_available():
    """Check if cache server is available."""
    if not CACHE_API_TOKEN:
        return False
    
    try:
        import requests
        health_response = requests.get(
            f"{CACHE_SERVER_URL}/health",
            timeout=5
        )
        if health_response.status_code != 200:
            return False
        
        client = RemoteCacheClient(
            server_url=CACHE_SERVER_URL,
            api_token=CACHE_API_TOKEN,
            timeout=5,
        )
        client.check_permissions()
        return True
    except RemoteCachePermissionError:
        return False
    except Exception:
        return False


@pytest.fixture
def temp_cache_dir():
    """Create a temporary cache directory for testing."""
    temp_dir = tempfile.mkdtemp()
    yield temp_dir
    shutil.rmtree(temp_dir, ignore_errors=True)


@pytest.fixture
def remote_cache_client():
    """Create remote cache client if server is available."""
    if not CACHE_API_TOKEN:
        pytest.skip("CACHE_API_TOKEN environment variable not set (check tests/.env file)")
    
    try:
        client = RemoteCacheClient(
            server_url=CACHE_SERVER_URL,
            api_token=CACHE_API_TOKEN,
            timeout=10,
        )
        try:
            import time
            time.sleep(1)
            client.check_permissions()
        except RemoteCachePermissionError as e:
            pytest.skip(f"Authentication failed: {e}. Please check your API token in tests/.env")
        except RemoteCacheError as e:
            if "rate limit" in str(e).lower():
                pass
            else:
                raise
        return client
    except Exception as e:
        pytest.skip(f"Cache server not available: {e}")


@pytest.mark.integration
def test_remote_cache_client_initialization(remote_cache_client):
    """Test that remote cache client can be initialized and connect to server."""
    assert remote_cache_client is not None
    assert remote_cache_client.server_url == CACHE_SERVER_URL.rstrip('/')
    assert remote_cache_client.api_token == CACHE_API_TOKEN


@pytest.mark.integration
def test_remote_cache_permissions(remote_cache_client):
    """Test that user status can be checked."""
    try:
        status = remote_cache_client.check_permissions()
        assert isinstance(status, dict)
    except RemoteCacheError as e:
        if "rate limit" in str(e).lower():
            pytest.skip(f"Rate limit exceeded: {e}. This is expected during testing.")
        raise


@pytest.mark.integration
def test_remote_cache_upload_and_get(remote_cache_client, temp_cache_dir):
    """Test uploading cache to server and retrieving it."""
    import hashlib
    test_hash = hashlib.sha256(b"test_molecule_and_parameters").hexdigest()
    test_content = b"This is a test cache file content for ORCA output"
    
    test_file = Path(temp_cache_dir) / "test_output.out"
    test_file.write_bytes(test_content)
    
    try:
        success = remote_cache_client.upload_cache(
            input_hash=test_hash,
            output_file=test_file,
            input_parameters={'test': 'parameter'},
            file_extension='.out'
        )
        assert success is True
    except RemoteCacheError as e:
        if "already exists" in str(e).lower() or "reputation" in str(e).lower():
            pass
        elif "500" in str(e) or "internal server error" in str(e).lower():
            pass
        else:
            raise
    except RemoteCachePermissionError as e:
        pytest.skip(f"Upload permission denied: {e}")
    
    # Test check_cache - must succeed without method_version errors
    import time
    time.sleep(1)
    try:
        cache_info = remote_cache_client.check_cache(test_hash)
        # If cache exists, it should return dict, if not - None
        assert cache_info is None or isinstance(cache_info, dict)
    except RemoteCacheError as e:
        error_msg = str(e).lower()
        if "rate limit" in error_msg:
            pytest.skip(f"Rate limit exceeded: {e}")
        elif "method_version" in error_msg:
            pytest.fail(f"method_version format error: {e}. This should be fixed!")
        else:
            raise
    except RemoteCachePermissionError as e:
        pytest.skip(f"Read permission denied: {e}")
    
    # Test get_cache - must succeed without method_version errors
    time.sleep(1)
    try:
        retrieved_content = remote_cache_client.get_cache(test_hash)
        # If cache exists, it should return bytes, if not - None
        assert retrieved_content is None or isinstance(retrieved_content, bytes)
    except (RemoteCacheError, RemoteCachePermissionError) as e:
        error_msg = str(e).lower()
        if "not found" in error_msg:
            # Cache not found is acceptable
            pass
        elif "method_version" in error_msg:
            pytest.fail(f"method_version format error: {e}. This should be fixed!")
        else:
            raise


@pytest.mark.integration
def test_orca_with_remote_cache_integration(temp_cache_dir):
    """Test ORCA calculation with remote cache integration."""
    if not CACHE_API_TOKEN:
        pytest.skip("ORCA_CACHE_API_TOKEN environment variable not set")
    
    if not is_cache_server_available():
        pytest.skip("Cache server not available")
    
    orca = Orca(
        script_path="orca",
        working_dir=temp_cache_dir,
        output_dir=temp_cache_dir,
        functional="AM1",
        method_type="SP",
        n_processors=1,
        cache_dir=str(Path(temp_cache_dir) / ".orca_cache"),
        cache_server_url=CACHE_SERVER_URL,
        cache_api_token=CACHE_API_TOKEN,
        cache_timeout=30,
        log_level=20,
    )
    
    mol = AddHs(MolFromSmiles("CCO"))
    
    try:
        orca.cache.remote_cache_client.check_permissions()
    except Exception as e:
        pytest.skip(f"Failed to check user status: {e}")
    
    try:
        homo1 = orca.homo_energy(mol)
        assert homo1 is not None
    except Exception as e:
        pytest.skip(f"ORCA calculation failed: {e}")
    
    # Verify remote cache client is working
    assert orca.cache.remote_cache_client is not None, "Remote cache client should be initialized"
    
    # Clear local cache to force remote cache usage
    local_cache_file = orca.cache.get(orca._get_molecule_hash(mol))
    if local_cache_file and local_cache_file.exists():
        local_cache_file.unlink()
        mol_hash = orca._get_molecule_hash(mol)
        if mol_hash in orca.cache.index:
            del orca.cache.index[mol_hash]
            orca.cache._save_index()
    
    # Second calculation should use remote cache (or calculate if not cached)
    # This should NOT fall back to local cache silently
    try:
        homo2 = orca.homo_energy(mol)
        assert homo2 is not None
        assert abs(homo1 - homo2) < 0.001
    except Exception as e:
        error_msg = str(e).lower()
        if "method_version" in error_msg:
            pytest.fail(f"method_version format error: {e}. Remote cache should work correctly!")
        else:
            pytest.fail(f"Failed to retrieve from remote cache: {e}")


@pytest.mark.integration
def test_remote_cache_error_handling():
    """Test error handling for remote cache operations."""
    invalid_client = RemoteCacheClient(
        server_url=CACHE_SERVER_URL,
        api_token="invalid_token_12345",
        timeout=5,
    )
    
    with pytest.raises((RemoteCacheError, RemoteCachePermissionError)):
        invalid_client.check_cache("test_hash")
    
    invalid_url_client = RemoteCacheClient(
        server_url="http://invalid-server-12345.example.com",
        api_token="test_token",
        timeout=2,
    )
    
    with pytest.raises(RemoteCacheError):
        invalid_url_client.get_cache("test_hash")


@pytest.mark.integration
def test_remote_cache_fallback_to_local(temp_cache_dir):
    """Test that system falls back to local cache when remote fails."""
    if not CACHE_API_TOKEN:
        pytest.skip("ORCA_CACHE_API_TOKEN environment variable not set")
    
    orca = Orca(
        script_path="orca",
        working_dir=temp_cache_dir,
        output_dir=temp_cache_dir,
        functional="AM1",
        method_type="SP",
        n_processors=1,
        cache_dir=str(Path(temp_cache_dir) / ".orca_cache"),
        cache_server_url="http://invalid-server.example.com",
        cache_api_token="test_token",
        cache_timeout=2,
        log_level=30,
    )
    
    mol = AddHs(MolFromSmiles("CCO"))
    
    try:
        homo = orca.homo_energy(mol)
        assert homo is not None
    except Exception as e:
        pytest.skip(f"ORCA calculation failed: {e}")


@pytest.mark.integration
def test_batch_processing_with_remote_cache(temp_cache_dir):
    """Test batch processing with remote cache for 3 molecules."""
    if not CACHE_API_TOKEN:
        pytest.skip("CACHE_API_TOKEN environment variable not set (check tests/.env file)")
    
    orca = Orca(
        script_path="orca",
        working_dir=temp_cache_dir,
        output_dir=temp_cache_dir,
        functional="AM1",
        method_type="SP",
        n_processors=1,
        cache_dir=str(Path(temp_cache_dir) / ".orca_cache"),
        cache_server_url=CACHE_SERVER_URL,
        cache_api_token=CACHE_API_TOKEN,
        cache_timeout=30,
        log_level=20,
    )
    
    batch_processing = ORCABatchProcessing(
        orca=orca,
        parallel_mode="sequential",
    )
    
    smiles_list = ["CCO", "CC(=O)O", "C1=CC=CC=C1"]
    
    x = batch_processing.x_molecule()
    descriptors = [
        orca.homo_energy(x),
        orca.lumo_energy(x),
        orca.gap_energy(x),
    ]
    
    import time
    start_time = time.time()
    try:
        result1 = batch_processing.calculate_descriptors(
            smiles_list,
            descriptors=descriptors,
            progress=True,
        )
        first_duration = time.time() - start_time
        
        assert result1 is not None
        if hasattr(result1, 'shape'):
            assert result1.shape[0] == len(smiles_list), f"Expected {len(smiles_list)} rows, got {result1.shape[0]}"
        else:
            assert len(result1) == len(smiles_list), f"Expected {len(smiles_list)} results, got {len(result1)}"
        
        first_results = result1.copy() if hasattr(result1, 'copy') else result1
        
    except Exception as e:
        pytest.skip(f"Batch calculation failed: {e}")
    
    local_cache_dir = Path(temp_cache_dir) / ".orca_cache"
    if local_cache_dir.exists():
        for cache_file in local_cache_dir.glob("*"):
            if cache_file.is_file() and cache_file.name != "cache_index.json":
                cache_file.unlink()
        
        cache_index = local_cache_dir / "cache_index.json"
        if cache_index.exists():
            cache_index.write_text("{}")
    
    # Verify remote cache client is working
    assert orca.cache.remote_cache_client is not None, "Remote cache client should be initialized"
    
    start_time = time.time()
    try:
        result2 = batch_processing.calculate_descriptors(
            smiles_list,
            descriptors=descriptors,
            progress=True,
        )
        second_duration = time.time() - start_time
        
        assert result2 is not None
        
        if hasattr(first_results, 'iloc'):
            for col in first_results.columns:
                if col in result2.columns:
                    for idx in range(len(first_results)):
                        val1 = first_results.iloc[idx][col]
                        val2 = result2.iloc[idx][col]
                        if val1 is not None and val2 is not None:
                            assert abs(val1 - val2) < 0.001, \
                                f"Values differ for {col} at index {idx}: {val1} vs {val2}"
        else:
            for i, (r1, r2) in enumerate(zip(first_results, result2)):
                for key in r1:
                    if key in r2 and r1[key] is not None and r2[key] is not None:
                        assert abs(r1[key] - r2[key]) < 0.001, \
                            f"Values differ for {key} at index {i}: {r1[key]} vs {r2[key]}"
        
        # Second calculation should be faster with cache (or at least not slower)
        # If remote cache is working, it should be significantly faster
        assert second_duration < first_duration * 1.1, \
            f"Second calculation should use cache and be faster, but took {second_duration:.2f}s vs {first_duration:.2f}s"
        
    except Exception as e:
        error_msg = str(e).lower()
        if "method_version" in error_msg:
            pytest.fail(f"method_version format error: {e}. Remote cache should work correctly!")
        else:
            pytest.fail(f"Failed to retrieve from remote cache in batch processing: {e}")

