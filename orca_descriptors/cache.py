"""Caching system for ORCA calculation results."""

import json
import shutil
from pathlib import Path
from typing import Optional


class CacheManager:
    """Manage cache for ORCA calculation results."""
    
    def __init__(self, cache_dir: str):
        """Initialize cache manager.
        
        Args:
            cache_dir: Directory for storing cached results
        """
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.index_file = self.cache_dir / "cache_index.json"
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
        
        Args:
            mol_hash: Hash of the molecule and calculation parameters
            
        Returns:
            Path to cached output file, or None if not found
        """
        if mol_hash in self.index:
            cached_path = Path(self.index[mol_hash])
            if cached_path.exists():
                return cached_path
            else:
                # Remove invalid entry
                del self.index[mol_hash]
                self._save_index()
        
        return None
    
    def store(self, mol_hash: str, output_file: Path):
        """Store output file in cache.
        
        Args:
            mol_hash: Hash of the molecule and calculation parameters
            output_file: Path to ORCA output file
        """
        # Copy file to cache directory, preserving original extension
        if output_file.exists():
            cached_file = self.cache_dir / f"{mol_hash}{output_file.suffix}"
            shutil.copy2(output_file, cached_file)
            self.index[mol_hash] = str(cached_file)
            self._save_index()
            return cached_file
        return None
    
    def clear(self):
        """Clear all cached files."""
        if self.cache_dir.exists():
            shutil.rmtree(self.cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.index = {}
        self._save_index()

