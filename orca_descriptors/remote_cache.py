"""Remote cache client for ORCA calculation results."""

import json
import logging
import re
import time
from pathlib import Path
from typing import Optional, Dict, Any
import requests
from requests.exceptions import RequestException

logger = logging.getLogger(__name__)


class RemoteCacheError(Exception):
    """Base exception for remote cache errors."""
    pass


class RemoteCachePermissionError(RemoteCacheError):
    """Exception raised when user doesn't have required permissions."""
    pass


class RemoteCacheClient:
    """Client for interacting with remote cache service API.
    
    Uses API v1 endpoints:
    - POST /api/v1/cache/check - Check if cache exists
    - POST /api/v1/cache/upload - Upload cache
    - GET /api/v1/cache/{cache_id} - Get cache info
    - GET /api/v1/cache/{cache_id}/files/{filename} - Download file
    """
    
    def __init__(
        self,
        server_url: str = "https://api.orca-descriptors.massonnn.ru",
        api_token: str = None,
        timeout: int = 30,
        method_version: str = "6.0.1",
        max_retries: int = 5,
        retry_backoff_factor: float = 2.0,
    ):
        """Initialize remote cache client.
        
        Args:
            server_url: Base URL of the cache service (default: "https://api.orca-descriptors.massonnn.ru")
            api_token: API token for authentication (required)
            timeout: Request timeout in seconds (default: 30)
            method_version: ORCA version string (default: "6.0.1", used only for upload, format: "X.Y.Z")
            max_retries: Maximum number of retries for rate-limited requests (default: 5)
            retry_backoff_factor: Backoff factor for exponential backoff (default: 2.0)
        """
        self.server_url = server_url.rstrip('/')
        self.api_token = api_token
        self.timeout = timeout
        self.max_retries = max_retries
        self.retry_backoff_factor = retry_backoff_factor
        self.session = requests.Session()
        if api_token:
            self.session.headers.update({
                'X-API-Key': api_token,
                'Authorization': f'Bearer {api_token}',
            })
        self.method_version = method_version
    
    @property
    def method_version(self) -> str:
        """Get normalized method version."""
        return self._method_version
    
    @method_version.setter
    def method_version(self, value: str) -> None:
        """Set method version with automatic normalization."""
        self._method_version = self._normalize_method_version(value or "6.0.1")
    
    def _normalize_method_version(self, version: str) -> str:
        """Normalize method_version to expected API format.
        
        API expects format like "5.0.4" or "6.0.1" (version number only, without "ORCA" prefix).
        Based on API documentation and testing, the format should be just the version number.
        
        Args:
            version: Version string (may be in various formats like "ORCA 5.0.4", "5.0.4", etc.)
            
        Returns:
            Normalized version string in format "X.Y.Z" (e.g., "6.0.1")
        """
        if not version:
            return "6.0.1"
        
        version = version.strip()
        
        match = re.search(r"(\d+\.\d+\.\d+)", version)
        if match:
            return match.group(1)
        
        return "6.0.1"
    
    def _make_request(
        self,
        method: str,
        endpoint: str,
        retry_on_rate_limit: bool = True,
        **kwargs
    ) -> requests.Response:
        """Make HTTP request to the API with automatic retry on rate limiting.
        
        Args:
            method: HTTP method (GET, POST, etc.)
            endpoint: API endpoint (e.g., "/api/v1/cache/check")
            retry_on_rate_limit: Whether to retry on rate limit (429) errors (default: True)
            **kwargs: Additional arguments for requests
            
        Returns:
            Response object
            
        Raises:
            RemoteCacheError: For general API errors
            RemoteCachePermissionError: For permission errors
        """
        url = f"{self.server_url}{endpoint}"
        
        for attempt in range(self.max_retries + 1):
            try:
                response = self.session.request(
                    method=method,
                    url=url,
                    timeout=self.timeout,
                    **kwargs
                )
                
                if response.status_code == 401:
                    raise RemoteCachePermissionError(
                        "Authentication failed. Please check your API token."
                    )
                elif response.status_code == 403:
                    raise RemoteCachePermissionError(
                        "Access denied. Please check your API token permissions."
                    )
                elif response.status_code == 404:
                    return response
                elif response.status_code == 400:
                    try:
                        error_data = response.json()
                        error_msg = error_data.get('detail', 'Bad Request')
                        if isinstance(error_msg, list):
                            error_msg = '; '.join([err.get('msg', str(err)) for err in error_msg])
                    except (ValueError, KeyError, json.JSONDecodeError):
                        error_msg = "Bad Request"
                    logger.debug(f"Bad Request error: {error_msg}. Skipping request.")
                    raise RemoteCacheError(f"Bad Request: {error_msg}")
                elif response.status_code == 429:
                    if retry_on_rate_limit and attempt < self.max_retries:
                        wait_time = self.retry_backoff_factor ** attempt
                        retry_after = response.headers.get('Retry-After')
                        if retry_after:
                            try:
                                wait_time = float(retry_after)
                            except (ValueError, TypeError):
                                pass
                        
                        logger.debug(
                            f"Rate limit exceeded (attempt {attempt + 1}/{self.max_retries + 1}). "
                            f"Waiting {wait_time:.2f} seconds before retry..."
                        )
                        time.sleep(wait_time)
                        continue
                    else:
                        try:
                            error_data = response.json()
                            if 'detail' in error_data:
                                if isinstance(error_data['detail'], list):
                                    error_msgs = [err.get('msg', str(err)) for err in error_data['detail']]
                                    error_msg = '; '.join(error_msgs)
                                else:
                                    error_msg = str(error_data['detail'])
                            else:
                                error_msg = error_data.get('error', error_data.get('message', 'Rate limit exceeded'))
                        except (ValueError, KeyError, json.JSONDecodeError):
                            error_msg = f"Rate limit exceeded: HTTP {response.status_code}"
                        raise RemoteCacheError(f"Rate limit exceeded after {attempt + 1} attempts: {error_msg}")
                elif response.status_code >= 400:
                    try:
                        error_data = response.json()
                        if 'detail' in error_data:
                            if isinstance(error_data['detail'], list):
                                error_msgs = [err.get('msg', str(err)) for err in error_data['detail']]
                                error_msg = '; '.join(error_msgs)
                            else:
                                error_msg = str(error_data['detail'])
                        else:
                            error_msg = error_data.get('error', error_data.get('message', 'Unknown error'))
                        logger.debug(f"API error response (status {response.status_code}): {json.dumps(error_data, indent=2)}")
                        raise RemoteCacheError(f"API error: {error_msg}")
                    except (ValueError, KeyError, json.JSONDecodeError):
                        error_text = response.text[:500]
                        logger.debug(f"API error response (status {response.status_code}): {error_text}")
                        raise RemoteCacheError(
                            f"API error: HTTP {response.status_code} - {error_text}"
                        )
                
                return response
                
            except requests.exceptions.Timeout:
                if attempt < self.max_retries:
                    wait_time = self.retry_backoff_factor ** attempt
                    logger.debug(f"Request timeout (attempt {attempt + 1}). Waiting {wait_time:.2f} seconds...")
                    time.sleep(wait_time)
                    continue
                raise RemoteCacheError(
                    f"Request timeout after {self.timeout} seconds and {attempt + 1} attempts. "
                    f"The server may be unavailable or too slow."
                )
            except requests.exceptions.ConnectionError as e:
                if attempt < self.max_retries:
                    wait_time = self.retry_backoff_factor ** attempt
                    logger.debug(f"Connection error (attempt {attempt + 1}). Waiting {wait_time:.2f} seconds...")
                    time.sleep(wait_time)
                    continue
                raise RemoteCacheError(
                    f"Failed to connect to cache server at {self.server_url} after {attempt + 1} attempts. "
                    f"Please check if the server is running and the URL is correct."
                )
            except RequestException as e:
                if attempt < self.max_retries:
                    wait_time = self.retry_backoff_factor ** attempt
                    logger.debug(f"Network error (attempt {attempt + 1}). Waiting {wait_time:.2f} seconds...")
                    time.sleep(wait_time)
                    continue
                raise RemoteCacheError(f"Network error after {attempt + 1} attempts: {str(e)}")
        
        raise RemoteCacheError(f"Request failed after {self.max_retries + 1} attempts")
    
    def check_cache(
        self,
        input_hash: str,
        only_validated: bool = False
    ) -> Optional[Dict[str, Any]]:
        """Check if cache exists for given input hash.
        
        Note: Cache is shared across all ORCA versions, so method_version is not used.
        
        Args:
            input_hash: Hash of the molecule and calculation parameters
            only_validated: If True, returns only validated caches (default: False)
            
        Returns:
            CacheResponse dict if cache found, None otherwise
            
        Raises:
            RemoteCacheError: For API errors
            RemoteCachePermissionError: If user doesn't have read permission
        """
        try:
            normalized_version = self._normalize_method_version(self.method_version or "6.0.1")
            request_payload = {
                'input_hash': input_hash,
                'method_version': normalized_version,
                'only_validated': only_validated,
            }
            logger.debug(f"Sending check_cache request with payload: {json.dumps(request_payload, indent=2)}")
            response = self._make_request(
                'POST',
                '/api/v1/cache/check',
                json=request_payload
            )
            
            if response.status_code == 404:
                return None
            
            response.raise_for_status()
            cache_data = response.json()
            
            if cache_data is None:
                return None
            
            return cache_data
            
        except RemoteCachePermissionError:
            raise
        except Exception as e:
            if isinstance(e, RemoteCacheError):
                raise
            raise RemoteCacheError(f"Failed to check cache: {str(e)}")
    
    def get_cache(self, input_hash: str) -> Optional[bytes]:
        """Retrieve cached output file from server.
        
        First checks if cache exists, then downloads the main output file.
        Uses file_paths from check_cache response to find the file.
        
        Args:
            input_hash: Hash of the molecule and calculation parameters
            
        Returns:
            Cached file content as bytes, or None if not found
            
        Raises:
            RemoteCacheError: For API errors
            RemoteCachePermissionError: If user doesn't have read permission
        """
        try:
            cache_info = self.check_cache(input_hash)
            if cache_info is None:
                return None
            
            cache_id = cache_info.get('id')
            if cache_id is None:
                return None
            
            file_paths = cache_info.get('file_paths', [])
            if not file_paths:
                logger.debug(f"Cache {cache_id} has no files")
                return None
            
            filename = None
            for file_path in file_paths:
                path_parts = file_path.split('/')
                file_name = path_parts[-1] if path_parts else file_path
                
                if file_name.endswith('.out') or file_name.endswith('.log') or file_name.endswith('.smd.out'):
                    filename = file_name
                    if file_name.endswith('.out'):
                        break
            
            if not filename:
                path_parts = file_paths[0].split('/')
                filename = path_parts[-1] if path_parts else file_paths[0]
            
            file_response = self._make_request(
                'GET',
                f'/api/v1/cache/{cache_id}/files/{filename}'
            )
            file_response.raise_for_status()
            
            return file_response.content
            
        except RemoteCachePermissionError:
            raise
        except Exception as e:
            if isinstance(e, RemoteCacheError):
                raise
            raise RemoteCacheError(f"Failed to retrieve cache: {str(e)}")
    
    def upload_cache(
        self,
        input_hash: str,
        output_file: Path,
        input_parameters: Optional[Dict[str, Any]] = None,
        file_extension: str = '.out'
    ) -> bool:
        """Upload output file to cache server.
        
        Args:
            input_hash: Hash of the molecule and calculation parameters
            output_file: Path to ORCA output file
            input_parameters: Dictionary with calculation parameters (optional)
            file_extension: File extension (default: '.out', not used but kept for compatibility)
            
        Returns:
            True if upload was successful, False otherwise
            
        Raises:
            RemoteCacheError: For API errors
            RemoteCachePermissionError: If user doesn't have upload permission
        """
        if not output_file.exists():
            logger.debug(f"Cannot upload non-existent file: {output_file}")
            return False
        
        if input_parameters is None:
            input_parameters = {}
        
        url = f"{self.server_url}/api/v1/cache/upload"
        
        for attempt in range(self.max_retries + 1):
            try:
                with open(output_file, 'rb') as f:
                    files = [
                        ('files', (output_file.name, f, 'application/octet-stream'))
                    ]
                    normalized_version = self._normalize_method_version(self.method_version or "6.0.1")
                    data = {
                        'input_hash': input_hash,
                        'method_version': normalized_version,
                        'input_parameters': json.dumps(input_parameters),
                    }
                    
                    headers = {k: v for k, v in self.session.headers.items() 
                              if k.lower() != 'content-type'}
                    
                    if 'X-API-Key' not in headers:
                        headers['X-API-Key'] = self.api_token
                    if 'Authorization' not in headers:
                        headers['Authorization'] = f'Bearer {self.api_token}'
                    
                    response = self.session.post(
                        url,
                        files=files,
                        data=data,
                        headers=headers,
                        timeout=self.timeout,
                    )
                    
                    if response.status_code == 401:
                        raise RemoteCachePermissionError(
                            "Authentication failed. Please check your API token."
                        )
                    elif response.status_code == 403:
                        raise RemoteCachePermissionError(
                            "Access denied. You don't have permission to upload cache. "
                            "Please check your API token permissions on the website."
                        )
                    elif response.status_code == 400:
                        try:
                            error_data = response.json()
                            error_msg = error_data.get('detail', 'Bad Request')
                            if isinstance(error_msg, list):
                                error_msg = '; '.join([err.get('msg', str(err)) for err in error_msg])
                        except (ValueError, KeyError, json.JSONDecodeError):
                            error_msg = "Bad Request"
                        logger.debug(f"Bad Request error during upload: {error_msg}. Skipping upload.")
                        raise RemoteCacheError(f"Bad Request: {error_msg}")
                    elif response.status_code == 429:
                        if attempt < self.max_retries:
                            wait_time = self.retry_backoff_factor ** attempt
                            retry_after = response.headers.get('Retry-After')
                            if retry_after:
                                try:
                                    wait_time = float(retry_after)
                                except (ValueError, TypeError):
                                    pass
                            logger.debug(
                                f"Rate limit exceeded during upload (attempt {attempt + 1}/{self.max_retries + 1}). "
                                f"Waiting {wait_time:.2f} seconds before retry..."
                            )
                            time.sleep(wait_time)
                            continue
                        else:
                            try:
                                error_data = response.json()
                                if 'detail' in error_data:
                                    if isinstance(error_data['detail'], list):
                                        error_msgs = [err.get('msg', str(err)) for err in error_data['detail']]
                                        error_msg = '; '.join(error_msgs)
                                    else:
                                        error_msg = str(error_data['detail'])
                                else:
                                    error_msg = error_data.get('error', error_data.get('message', 'Rate limit exceeded'))
                            except (ValueError, KeyError, json.JSONDecodeError):
                                error_msg = f"Rate limit exceeded: HTTP {response.status_code}"
                            raise RemoteCacheError(f"Rate limit exceeded after {attempt + 1} attempts: {error_msg}")
                    elif response.status_code == 201:
                        return True
                    elif response.status_code >= 400:
                        try:
                            error_data = response.json()
                            if 'detail' in error_data:
                                if isinstance(error_data['detail'], list):
                                    error_msgs = [err.get('msg', str(err)) for err in error_data['detail']]
                                    error_msg = '; '.join(error_msgs)
                                else:
                                    error_msg = str(error_data['detail'])
                            else:
                                error_msg = error_data.get('error', error_data.get('message', 'Unknown error'))
                            raise RemoteCacheError(f"API error: {error_msg}")
                        except (ValueError, KeyError, json.JSONDecodeError):
                            raise RemoteCacheError(
                                f"API error: HTTP {response.status_code} - {response.text[:200]}"
                            )
                    
                    response.raise_for_status()
                    return True
                    
            except RemoteCachePermissionError:
                raise
            except requests.exceptions.Timeout:
                if attempt < self.max_retries:
                    wait_time = self.retry_backoff_factor ** attempt
                    logger.debug(f"Upload timeout (attempt {attempt + 1}). Waiting {wait_time:.2f} seconds...")
                    time.sleep(wait_time)
                    continue
                raise RemoteCacheError(
                    f"Upload timeout after {self.timeout} seconds and {attempt + 1} attempts. "
                    f"The server may be unavailable or too slow."
                )
            except requests.exceptions.ConnectionError as e:
                if attempt < self.max_retries:
                    wait_time = self.retry_backoff_factor ** attempt
                    logger.debug(f"Connection error during upload (attempt {attempt + 1}). Waiting {wait_time:.2f} seconds...")
                    time.sleep(wait_time)
                    continue
                raise RemoteCacheError(
                    f"Failed to connect to cache server at {self.server_url} after {attempt + 1} attempts. "
                    f"Please check if the server is running and the URL is correct."
                )
            except RequestException as e:
                if attempt < self.max_retries:
                    wait_time = self.retry_backoff_factor ** attempt
                    logger.debug(f"Network error during upload (attempt {attempt + 1}). Waiting {wait_time:.2f} seconds...")
                    time.sleep(wait_time)
                    continue
                raise RemoteCacheError(f"Network error during upload after {attempt + 1} attempts: {str(e)}")
            except Exception as e:
                if isinstance(e, (RemoteCacheError, RemoteCachePermissionError)):
                    raise
                if attempt < self.max_retries:
                    wait_time = self.retry_backoff_factor ** attempt
                    logger.debug(f"Error during upload (attempt {attempt + 1}). Waiting {wait_time:.2f} seconds...")
                    time.sleep(wait_time)
                    continue
                raise RemoteCacheError(f"Failed to upload cache after {attempt + 1} attempts: {str(e)}")
        
        raise RemoteCacheError(f"Upload failed after {self.max_retries + 1} attempts")
    
    def check_permissions(self) -> dict:
        """Check user status and permissions.
        
        Note: This endpoint may not be available in all API versions.
        If the endpoint is not found, returns an empty dict.
        
        Returns:
            Dictionary with user status information, or empty dict if endpoint not available
            
        Raises:
            RemoteCacheError: For API errors (except 404)
            RemoteCachePermissionError: For permission errors
        """
        try:
            response = self._make_request('GET', '/api/v1/user/status')
            if response.status_code == 404:
                logger.debug("User status endpoint not available")
                return {}
            response.raise_for_status()
            return response.json()
        except RemoteCachePermissionError:
            raise
        except Exception as e:
            if isinstance(e, RemoteCacheError):
                raise
            if hasattr(e, 'response') and hasattr(e.response, 'status_code') and e.response.status_code == 404:
                return {}
            raise RemoteCacheError(f"Failed to check user status: {str(e)}")
