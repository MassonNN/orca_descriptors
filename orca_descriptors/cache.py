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
                    data = json.load(f)
                    # Handle both old format (dict of hash -> path) and new format (dict of hash -> dict)
                    if isinstance(data, dict) and data:
                        # Check if it's old format (all values are strings) or new format
                        first_value = next(iter(data.values()))
                        if isinstance(first_value, str):
                            # Old format: convert to new format
                            self.index = {k: {'path': v, 'input_parameters': None} for k, v in data.items()}
                        else:
                            # New format
                            self.index = data
                    else:
                        self.index = {}
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
        if mol_hash in self.index:
            index_entry = self.index[mol_hash]
            if isinstance(index_entry, str):
                cached_path = Path(index_entry)
                input_parameters = None
            else:
                cached_path = Path(index_entry.get('path', ''))
                input_parameters = index_entry.get('input_parameters')
            
            if cached_path.exists():
                # Try to extract input_parameters from output file if not available in index
                if not input_parameters or not isinstance(input_parameters, dict) or len(input_parameters) == 0:
                    try:
                        from orca_descriptors.output_parser import ORCAOutputParser
                        parser = ORCAOutputParser()
                        parsed_data = parser.parse(cached_path)
                        extracted_params = parsed_data.get('input_parameters', {})
                        if extracted_params and isinstance(extracted_params, dict) and len(extracted_params) > 0:
                            input_parameters = extracted_params
                            logger.debug(f"Extracted input_parameters from output file for hash: {mol_hash}")
                    except Exception as e:
                        logger.debug(f"Could not extract input_parameters from output file: {e}")
                
                # Try to extract from input file if still not available
                if not input_parameters or not isinstance(input_parameters, dict) or len(input_parameters) == 0:
                    try:
                        # Try to find corresponding input file
                        input_file = cached_path.parent / cached_path.name.replace('.out', '.inp').replace('.log', '.inp').replace('.smd.out', '.inp')
                        if not input_file.exists():
                            # Try in working directory
                            base_name = cached_path.stem.replace('.smd', '')
                            input_file = cached_path.parent.parent / f"{base_name}.inp"
                        
                        if input_file.exists():
                            from orca_descriptors.output_parser import ORCAOutputParser
                            parser = ORCAOutputParser()
                            extracted_params = parser.parse_input_file(input_file)
                            if extracted_params and isinstance(extracted_params, dict) and len(extracted_params) > 0:
                                input_parameters = extracted_params
                                logger.debug(f"Extracted input_parameters from input file for hash: {mol_hash}")
                    except Exception as e:
                        logger.debug(f"Could not extract input_parameters from input file: {e}")
                
                # Update index with extracted parameters if found
                if input_parameters and isinstance(input_parameters, dict) and len(input_parameters) > 0:
                    if isinstance(index_entry, str):
                        self.index[mol_hash] = {
                            'path': str(cached_path),
                            'input_parameters': input_parameters
                        }
                    else:
                        self.index[mol_hash]['input_parameters'] = input_parameters
                    self._save_index()
                
                if self.remote_cache_client:
                    try:
                        logger.debug(f"Found local cache for hash: {mol_hash}, checking remote cache...")
                        orca_version = None
                        try:
                            from orca_descriptors.output_parser import ORCAOutputParser
                            parser = ORCAOutputParser()
                            parsed_data = parser.parse(cached_path)
                            orca_version = parsed_data.get('orca_version')
                        except Exception as e:
                            logger.debug(f"Could not parse ORCA version from local cache: {e}")
                        
                        if orca_version:
                            self.remote_cache_client.method_version = self.remote_cache_client._normalize_method_version(orca_version)
                            logger.debug(f"Using ORCA version from local cache: {self.remote_cache_client.method_version}")
                        
                        cache_exists = False
                        try:
                            existing_cache = self.remote_cache_client.check_cache(mol_hash)
                            if existing_cache is not None:
                                cache_exists = True
                                logger.debug(f"Cache already exists on server for hash: {mol_hash}")
                        except Exception as check_error:
                            logger.debug(f"Could not check if cache exists on server: {check_error}")
                        
                        # Try to upload to remote cache if it doesn't exist (don't fail if it doesn't work)
                        if not cache_exists:
                            # Only upload if we have input_parameters (skip if None or empty)
                            if input_parameters and isinstance(input_parameters, dict) and len(input_parameters) > 0:
                                try:
                                    logger.debug(f"Uploading local cache to remote for hash: {mol_hash}")
                                    self.remote_cache_client.upload_cache(
                                        input_hash=mol_hash,
                                        output_file=cached_path,
                                        input_parameters=input_parameters,
                                        file_extension=cached_path.suffix
                                    )
                                    logger.debug(f"Successfully uploaded local cache to remote for hash: {mol_hash}")
                                except Exception as upload_error:
                                    error_msg = str(upload_error).lower()
                                    if "already exists" in error_msg or "reputation" in error_msg:
                                        logger.debug(f"Cache already exists on server or reputation issue: {upload_error}")
                                    elif "bad request" in error_msg:
                                        logger.debug(
                                            f"Bad Request error during upload: {upload_error}. "
                                            f"Skipping upload. Local cache is still available."
                                        )
                                    else:
                                        logger.debug(f"Failed to upload local cache to remote for hash {mol_hash}: {upload_error}")
                            else:
                                logger.debug(f"Skipping upload for hash {mol_hash}: input_parameters not available or empty")
                        else:
                            logger.debug(f"Skipping upload - cache already exists on server for hash: {mol_hash}")
                    except Exception as e:
                        # Log but don't fail - local cache is still available
                        logger.debug(f"Error while trying to upload local cache to remote for hash {mol_hash}: {e}")
                
                return cached_path
            else:
                # Remove invalid entry
                del self.index[mol_hash]
                self._save_index()
        
        # Check remote cache BEFORE local cache miss (to ensure we check both)
        # This is already done above, but we ensure it's checked even if local cache doesn't exist
        
        # Check remote cache if available
        if self.remote_cache_client:
            try:
                logger.debug(f"Checking remote cache for hash: {mol_hash}")
                # First check cache to get metadata including input_parameters
                cache_info = self.remote_cache_client.check_cache(mol_hash)
                
                if cache_info is not None:
                    # Extract input_parameters from cache_info
                    input_parameters = cache_info.get('input_parameters')
                    if input_parameters is None:
                        input_parameters = {}
                    
                    # Get file content
                    remote_content = self.remote_cache_client.get_cache(mol_hash)
                    
                    if remote_content is not None:
                        # Determine file extension from content or use default
                        # Try to detect from common ORCA output extensions
                        file_extension = '.out'
                        file_paths = cache_info.get('file_paths', [])
                        if file_paths:
                            # Try to get extension from file_paths
                            for file_path in file_paths:
                                path_parts = file_path.split('/')
                                file_name = path_parts[-1] if path_parts else file_path
                                if file_name.endswith('.out') or file_name.endswith('.log') or file_name.endswith('.smd.out'):
                                    file_extension = Path(file_name).suffix
                                    if file_name.endswith('.out'):
                                        break
                        else:
                            # Fallback: check if we have a record of the extension in local index
                            for ext in ['.out', '.log', '.smd.out']:
                                if mol_hash in self.index:
                                    index_entry = self.index[mol_hash]
                                    if isinstance(index_entry, str):
                                        old_path = Path(index_entry)
                                    else:
                                        old_path = Path(index_entry.get('path', ''))
                                    if old_path.suffix:
                                        file_extension = old_path.suffix
                                        break
                        
                        # Save to local cache
                        cached_file = self.cache_dir / f"{mol_hash}{file_extension}"
                        cached_file.write_bytes(remote_content)
                        # Store path and input_parameters in index (from remote cache)
                        self.index[mol_hash] = {
                            'path': str(cached_file),
                            'input_parameters': input_parameters  # Retrieved from remote cache
                        }
                        self._save_index()
                    
                        # Update method_version in remote cache client from cache_info or downloaded file
                        method_version = cache_info.get('method_version')
                        if method_version:
                            # Normalize version format before setting
                            self.remote_cache_client.method_version = self.remote_cache_client._normalize_method_version(method_version)
                            logger.debug(f"Updated method_version from cache_info: {self.remote_cache_client.method_version}")
                        else:
                            # Fallback: try to parse from downloaded file
                            try:
                                from orca_descriptors.output_parser import ORCAOutputParser
                                parser = ORCAOutputParser()
                                parsed_data = parser.parse(cached_file)
                                orca_version = parsed_data.get('orca_version')
                                if orca_version:
                                    # Normalize version format before setting
                                    self.remote_cache_client.method_version = self.remote_cache_client._normalize_method_version(orca_version)
                                    logger.debug(f"Updated method_version from downloaded cache: {self.remote_cache_client.method_version}")
                            except Exception as e:
                                logger.debug(f"Could not parse ORCA version from downloaded cache: {e}")
                        
                        logger.debug(f"Downloaded cache from remote: {cached_file} (with input_parameters: {input_parameters})")
                        return cached_file
                    
            except Exception as e:
                error_msg = str(e).lower()
                if "bad request" in error_msg:
                    logger.debug(
                        f"Bad Request error from remote cache: {e}. "
                        f"Skipping remote cache. Falling back to local cache only."
                    )
                else:
                    logger.debug(
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
            # Store both path and input_parameters in index
            self.index[mol_hash] = {
                'path': str(cached_file),
                'input_parameters': input_parameters if input_parameters is not None else {}
            }
            self._save_index()
            
            # Upload to remote cache if available
            if self.remote_cache_client:
                try:
                    logger.debug(f"Storing cache and uploading to remote for hash: {mol_hash}")
                    # Parse ORCA version from output file
                    orca_version = None
                    try:
                        from orca_descriptors.output_parser import ORCAOutputParser
                        parser = ORCAOutputParser()
                        parsed_data = parser.parse(cached_file)
                        orca_version = parsed_data.get('orca_version')
                    except Exception as e:
                        logger.debug(f"Could not parse ORCA version: {e}")
                    
                    # Update method_version in remote cache client if version was found
                    if orca_version:
                        # Normalize version format before setting
                        self.remote_cache_client.method_version = self.remote_cache_client._normalize_method_version(orca_version)
                        logger.debug(f"Using ORCA version from output: {self.remote_cache_client.method_version}")
                    
                    # Check if cache already exists on server before uploading
                    cache_exists = False
                    try:
                        existing_cache = self.remote_cache_client.check_cache(mol_hash)
                        if existing_cache is not None:
                            cache_exists = True
                            logger.debug(f"Cache already exists on server for hash: {mol_hash}, skipping upload")
                    except Exception as check_error:
                        # If check fails, try to upload anyway
                        logger.debug(f"Could not check if cache exists on server: {check_error}")
                    
                    # Always upload after calculation to ensure cache is on server
                    # Even if cache exists, we try to upload (server will reject if reputation is lower)
                    try:
                        logger.debug(f"Uploading cache to remote for hash: {mol_hash} (after calculation)")
                        # Ensure input_parameters is always a dict, never None
                        upload_params = input_parameters if input_parameters is not None else {}
                        self.remote_cache_client.upload_cache(
                            input_hash=mol_hash,
                            output_file=cached_file,
                            input_parameters=upload_params,
                            file_extension=output_file.suffix
                        )
                        logger.debug(f"Successfully uploaded cache to remote for hash: {mol_hash}")
                    except Exception as upload_error:
                        error_msg = str(upload_error).lower()
                        # Don't log as error if cache already exists or reputation issue (expected behavior)
                        if "already exists" in error_msg or "reputation" in error_msg:
                            logger.debug(f"Cache already exists on server or reputation issue: {upload_error}")
                        elif "bad request" in error_msg:
                            logger.debug(
                                f"Bad Request error during upload: {upload_error}. "
                                f"Skipping upload. Local cache is still available."
                            )
                        else:
                            logger.debug(
                                f"Failed to upload to remote cache for hash {mol_hash}: {upload_error}. "
                                f"Local cache is still available."
                            )
                except Exception as e:
                    # Log error but don't fail - local cache is still available
                    logger.debug(
                        f"Failed to upload to remote cache for hash {mol_hash}: {e}. "
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

