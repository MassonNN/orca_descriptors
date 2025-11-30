// Language switcher for ORCA Descriptors documentation

(function() {
    'use strict';
    
    var switcherCreated = false;
    var observer = null;
    
    function getCurrentLanguage() {
        var path = window.location.pathname;
        if (path.includes('/ru/')) {
            return 'ru';
        }
        if (path.includes('/en/')) {
            return 'en';
        }
        // Default to English if path doesn't contain language
        return 'en';
    }
    
    function getOtherLanguage() {
        return getCurrentLanguage() === 'ru' ? 'en' : 'ru';
    }
    
    function switchLanguage(e) {
        if (e) {
            e.preventDefault();
        }
        
        var currentLang = getCurrentLanguage();
        var otherLang = getOtherLanguage();
        var currentPath = window.location.pathname;
        var currentHash = window.location.hash;
        var currentSearch = window.location.search;
        
        // Normalize path (remove trailing slash except for root)
        var normalizedPath = currentPath;
        if (normalizedPath !== '/' && normalizedPath.endsWith('/')) {
            normalizedPath = normalizedPath.slice(0, -1);
        }
        
        var newPath = normalizedPath;
        
        // Handle language replacement
        // Structure: English at root (/), Russian at /ru/
        if (normalizedPath.startsWith('/ru/')) {
            // Currently on Russian version
            if (otherLang === 'en') {
                // Switch to English: remove /ru/ prefix
                newPath = normalizedPath.replace('/ru/', '/');
                // Handle case when /ru/ becomes empty
                if (newPath === '' || newPath === '/ru') {
                    newPath = '/';
                }
            } else {
                // Should not happen, but handle it
                newPath = normalizedPath.replace('/ru/', '/' + otherLang + '/');
            }
        } else if (normalizedPath === '/ru') {
            // Special case: /ru without trailing slash
            newPath = otherLang === 'en' ? '/' : '/ru/';
        } else if (normalizedPath.startsWith('/en/')) {
            // Currently on /en/ path (shouldn't happen in production, but handle it)
            if (otherLang === 'ru') {
                newPath = normalizedPath.replace('/en/', '/ru/');
            } else {
                newPath = normalizedPath.replace('/en/', '/');
            }
        } else if (normalizedPath === '/' || normalizedPath === '') {
            // Currently on root (English)
            if (otherLang === 'ru') {
                newPath = '/ru/';
            } else {
                newPath = '/';
            }
        } else {
            // Path doesn't have language prefix
            var pathParts = normalizedPath.split('/').filter(function(p) { return p; });
            var langIndex = -1;
            
            // Look for 'en' or 'ru' in path parts
            for (var i = 0; i < pathParts.length; i++) {
                if (pathParts[i] === 'en' || pathParts[i] === 'ru') {
                    langIndex = i;
                    break;
                }
            }
            
            if (langIndex >= 0) {
                // Found language in path, replace it
                pathParts[langIndex] = otherLang;
                newPath = '/' + pathParts.join('/');
            } else {
                // No language in path - this is English version (at root)
                if (otherLang === 'ru') {
                    // Switch to Russian: add /ru/ prefix
                    newPath = '/ru' + normalizedPath;
                } else {
                    // Stay at root
                    newPath = normalizedPath;
                }
            }
        }
        
        // Ensure path starts with /
        if (!newPath.startsWith('/')) {
            newPath = '/' + newPath;
        }
        
        // Preserve hash and search params
        var newUrl = newPath + currentSearch + currentHash;
        window.location.href = newUrl;
    }
    
    function removeExistingSwitchers() {
        var existing = document.querySelector('#language-switcher');
        if (existing) {
            existing.remove();
        }
        var existingTop = document.querySelector('#language-switcher-top');
        if (existingTop) {
            existingTop.remove();
        }
    }
    
    function createLanguageSwitcher() {
        // Prevent duplicate creation
        if (switcherCreated && document.querySelector('#language-switcher')) {
            return;
        }
        
        // Remove any existing switchers first
        removeExistingSwitchers();
        
        var currentLang = getCurrentLanguage();
        var otherLang = getOtherLanguage();
        var langNames = {
            'en': 'English',
            'ru': 'Русский'
        };
        
        // Try to find sidebar for Furo theme
        var sidebar = document.querySelector('.sidebar-container') || 
                     document.querySelector('aside.sidebar') ||
                     document.querySelector('aside') ||
                     document.querySelector('.sidebar');
        
        if (sidebar) {
            // Create language switcher button
            var switcher = document.createElement('div');
            switcher.id = 'language-switcher';
            
            var button = document.createElement('button');
            button.textContent = langNames[otherLang];
            button.className = 'language-switcher-btn';
            button.onclick = switchLanguage;
            button.setAttribute('type', 'button');
            
            switcher.appendChild(button);
            
            // Find the best place to insert in Furo sidebar
            var sidebarBrand = sidebar.querySelector('.sidebar-brand');
            var sidebarSearch = sidebar.querySelector('.sidebar-search');
            var sidebarScroll = sidebar.querySelector('.sidebar-scroll');
            
            if (sidebarBrand) {
                // Insert after brand
                if (sidebarBrand.nextSibling) {
                    sidebarBrand.parentNode.insertBefore(switcher, sidebarBrand.nextSibling);
                } else {
                    sidebarBrand.parentNode.appendChild(switcher);
                }
            } else if (sidebarSearch) {
                // Insert after search
                if (sidebarSearch.nextSibling) {
                    sidebarSearch.parentNode.insertBefore(switcher, sidebarSearch.nextSibling);
                } else {
                    sidebarSearch.parentNode.appendChild(switcher);
                }
            } else if (sidebarScroll) {
                // Insert at the beginning of scroll container
                sidebarScroll.insertBefore(switcher, sidebarScroll.firstChild);
            } else {
                // Insert at the beginning of sidebar
                sidebar.insertBefore(switcher, sidebar.firstChild);
            }
            
            switcherCreated = true;
        }
        
        // Also try to add to top navigation bar if exists
        var topNav = document.querySelector('.topbar') || 
                    document.querySelector('nav') ||
                    document.querySelector('header');
        if (topNav && !document.querySelector('#language-switcher-top')) {
            var topSwitcher = document.createElement('div');
            topSwitcher.id = 'language-switcher-top';
            topSwitcher.style.cssText = 'margin-left: auto; margin-right: 10px;';
            
            var topButton = document.createElement('a');
            topButton.href = '#';
            topButton.textContent = '🌐 ' + langNames[otherLang];
            topButton.className = 'language-switcher-link';
            topButton.onclick = function(e) {
                switchLanguage(e);
            };
            
            topSwitcher.appendChild(topButton);
            topNav.appendChild(topSwitcher);
        }
    }
    
    function initializeSwitcher() {
        // Wait for DOM to be ready
        if (document.readyState === 'loading') {
            document.addEventListener('DOMContentLoaded', initializeSwitcher);
            return;
        }
        
        // Wait for theme to initialize
        setTimeout(function() {
            createLanguageSwitcher();
            
            // Watch for theme changes (Furo theme switcher)
            if (!observer) {
                observer = new MutationObserver(function(mutations) {
                    var switcherExists = document.querySelector('#language-switcher');
                    if (!switcherExists || !switcherExists.parentNode) {
                        switcherCreated = false;
                        setTimeout(createLanguageSwitcher, 100);
                    }
                });
                
                // Observe sidebar for changes
                var sidebar = document.querySelector('.sidebar-container') || 
                             document.querySelector('aside.sidebar') ||
                             document.querySelector('aside');
                if (sidebar) {
                    observer.observe(sidebar, {
                        childList: true,
                        subtree: true
                    });
                }
                
                // Also observe body for theme changes
                observer.observe(document.body, {
                    attributes: true,
                    attributeFilter: ['class', 'data-theme']
                });
            }
        }, 300);
    }
    
    // Listen for theme changes
    document.addEventListener('DOMContentLoaded', function() {
        // Furo theme switcher may trigger events
        var themeButton = document.querySelector('button[data-theme-toggle]') ||
                         document.querySelector('[data-theme-toggle]');
        if (themeButton) {
            themeButton.addEventListener('click', function() {
                setTimeout(function() {
                    switcherCreated = false;
                    createLanguageSwitcher();
                }, 200);
            });
        }
    });
    
    // Initialize
    initializeSwitcher();
    
    // Re-initialize on page visibility change (handles theme persistence)
    document.addEventListener('visibilitychange', function() {
        if (!document.hidden) {
            setTimeout(function() {
                if (!document.querySelector('#language-switcher')) {
                    switcherCreated = false;
                    createLanguageSwitcher();
                }
            }, 100);
        }
    });
})();

