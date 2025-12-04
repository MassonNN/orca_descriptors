"""Caching system for ORCA calculation results."""

import json
import logging
import shutil
from pathlib import Path
from typing import Optional, Dict, Any

logger = logging.getLogger(__name__)


class CacheManager:
    """Manage cache for ORCA calculation results.
    
    Supports both local and remote caching. If remote_cache_client is provided,
    the cache manager will check remote cache if local cache misses, and upload
    to remote cache after storing locally.
    """
    
    def __init__(
        self,
        cache_dir: str,
        remote_cache_client: Optional[object] = None,
    ):
        """Initialize cache manager.
        
        Args:
            cache_dir: Directory for storing cached results
            remote_cache_client: Optional RemoteCacheClient instance for remote caching
        """
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.index_file = self.cache_dir / "cache_index.json"
        self.remote_cache_client = remote_cache_client
        self._load_index()
    
    def _load_index(self):
        """Load cache index from disk."""
        if self.index_file.exists():
            try:
                with open(self.index_file, "r") as f:
                    self.index = json.load(f)
            except (json.JSONDecodeError, IOError):
                self.index = {}
        else:
            self.index = {}
    
    def _save_index(self):
        """Save cache index to disk."""
        with open(self.index_file, "w") as f:
            json.dump(self.index, f, indent=2)
    
    def get(self, mol_hash: str) -> Optional[Path]:
        """Get cached output file path if it exists.
        
        Checks local cache first, then remote cache if available.
        If found in remote cache, downloads and stores locally.
        
        Args:
            mol_hash: Hash of the molecule and calculation parameters
            
        Returns:
            Path to cached output file, or None if not found
        """
        # Check local cache first
        if mol_hash in self.index:
            cached_path = Path(self.index[mol_hash])
            if cached_path.exists():
                return cached_path
            else:
                # Remove invalid entry
                del self.index[mol_hash]
                self._save_index()
        
        # Check remote cache if available
        if self.remote_cache_client:
            try:
                logger.debug(f"Checking remote cache for hash: {mol_hash}")
                # Use mol_hash as input_hash for API
                remote_content = self.remote_cache_client.get_cache(mol_hash)
                
                if remote_content is not None:
                    # Determine file extension from content or use default
                    # Try to detect from common ORCA output extensions
                    file_extension = '.out'
                    for ext in ['.out', '.log', '.smd.out']:
                        if mol_hash in self.index:
                            # Check if we have a record of the extension
                            old_path = Path(self.index[mol_hash])
                            if old_path.suffix:
                                file_extension = old_path.suffix
                                break
                    
                    # Save to local cache
                    cached_file = self.cache_dir / f"{mol_hash}{file_extension}"
                    cached_file.write_bytes(remote_content)
                    self.index[mol_hash] = str(cached_file)
                    self._save_index()
                    
                    logger.debug(f"Downloaded cache from remote: {cached_file}")
                    return cached_file
                    
            except Exception as e:
                # Log error but don't fail - fall back to local-only behavior
                logger.warning(
                    f"Failed to retrieve from remote cache: {e}. "
                    f"Continuing with local cache only."
                )
        
        return None
    
    def store(
        self,
        mol_hash: str,
        output_file: Path,
        input_parameters: Optional[dict] = None
    ):
        """Store output file in cache.
        
        Stores locally first, then uploads to remote cache if available.
        
        Args:
            mol_hash: Hash of the molecule and calculation parameters (used as input_hash)
            output_file: Path to ORCA output file
            input_parameters: Optional dictionary with calculation parameters for remote cache
        """
        # Copy file to cache directory, preserving original extension
        if output_file.exists():
            cached_file = self.cache_dir / f"{mol_hash}{output_file.suffix}"
            shutil.copy2(output_file, cached_file)
            self.index[mol_hash] = str(cached_file)
            self._save_index()
            
            # Upload to remote cache if available
            if self.remote_cache_client:
                try:
                    logger.debug(f"Uploading cache to remote for hash: {mol_hash}")
                    self.remote_cache_client.upload_cache(
                        input_hash=mol_hash,
                        output_file=cached_file,
                        input_parameters=input_parameters,
                        file_extension=output_file.suffix
                    )
                    logger.debug(f"Successfully uploaded cache to remote")
                except Exception as e:
                    # Log error but don't fail - local cache is still available
                    logger.warning(
                        f"Failed to upload to remote cache: {e}. "
                        f"Local cache is still available."
                    )
            
            return cached_file
        return None
    
    def clear(self):
        """Clear all cached files."""
        if self.cache_dir.exists():
            shutil.rmtree(self.cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.index = {}
        self._save_index()

