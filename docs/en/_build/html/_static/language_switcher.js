// Language switcher for ORCA Descriptors documentation

(function() {
    'use strict';
    
    var switcherCreated = false;
    var observer = null;
    
    function getCurrentLanguage() {
        // Check document language attribute
        var htmlLang = document.documentElement.getAttribute('lang');
        if (htmlLang) {
            var lang = htmlLang.split('-')[0].toLowerCase();
            if (lang === 'ru' || lang === 'en') {
                return lang;
            }
        }
        
        // Check path for language indicator
        var path = window.location.pathname;
        var href = window.location.href;
        
        // Check for /ru/ in path (production structure: English at root, Russian at /ru/)
        if (path.includes('/ru/') || path.endsWith('/ru') || href.includes('/ru/')) {
            return 'ru';
        }
        
        // Check for /en/ in path (if explicitly set)
        if (path.includes('/en/') || path.endsWith('/en') || href.includes('/en/')) {
            return 'en';
        }
        
        // Default to English (root path)
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
        var currentProtocol = window.location.protocol;
        var currentHref = window.location.href;
        var currentPath = window.location.pathname;
        var currentHash = window.location.hash;
        var currentSearch = window.location.search;
        
        var newUrl;
        
        if (currentProtocol === 'file:') {
            // For local files (file://), work with the full href
            // Structure: English at /en/_build/html/..., Russian at /ru/_build/html/...
            if (currentHref.includes('/en/_build/html/')) {
                newUrl = currentHref.replace('/en/_build/html/', '/ru/_build/html/');
            } else if (currentHref.includes('/ru/_build/html/')) {
                newUrl = currentHref.replace('/ru/_build/html/', '/en/_build/html/');
            } else if (currentHref.includes('/en/')) {
                newUrl = currentHref.replace('/en/', '/ru/');
            } else if (currentHref.includes('/ru/')) {
                newUrl = currentHref.replace('/ru/', '/en/');
            } else {
                // Fallback: add /ru/ or /en/ based on target language
                var baseUrl = currentHref.substring(0, currentHref.lastIndexOf('/'));
                var fileName = currentHref.substring(currentHref.lastIndexOf('/'));
                if (otherLang === 'ru') {
                    newUrl = baseUrl.replace(/\/en(\/|$)/, '/ru/') + fileName;
                } else {
                    newUrl = baseUrl.replace(/\/ru(\/|$)/, '/en/') + fileName;
                }
            }
        } else {
            // For http/https, work with pathname
            var currentHost = window.location.host;
            var normalizedPath = currentPath;
            if (normalizedPath !== '/' && normalizedPath.endsWith('/')) {
                normalizedPath = normalizedPath.slice(0, -1);
            }
            
            var newPath = normalizedPath;
            
            // Handle language replacement
            // Structure: /... (English) and /ru/... (Russian) in production
            if (normalizedPath.startsWith('/ru/')) {
                // Currently on Russian version (/ru/...)
                // Switch to English: remove /ru/ prefix
                newPath = normalizedPath.replace(/^\/ru\//, '/');
                if (newPath === '') {
                    newPath = '/';
                }
            } else if (normalizedPath === '/ru') {
                // Special case: /ru without trailing slash
                newPath = '/';
            } else {
                // Currently on English version (root)
                // Switch to Russian: add /ru/ prefix
                if (normalizedPath === '/') {
                    newPath = '/ru/';
                } else {
                    newPath = '/ru' + normalizedPath;
                }
            }
            
            // Ensure path starts with /
            if (!newPath.startsWith('/')) {
                newPath = '/' + newPath;
            }
            
            // Build URL for http/https
            newUrl = currentProtocol + '//' + currentHost + newPath + currentSearch + currentHash;
        }
        
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
            
            // Watch for sidebar structure changes
            if (!observer) {
                var sidebarObserver = new MutationObserver(function(mutations) {
                    var switcherExists = document.querySelector('#language-switcher');
                    if (!switcherExists || !switcherExists.parentNode) {
                        switcherCreated = false;
                        setTimeout(createLanguageSwitcher, 100);
                    }
                });
                
                // Observe sidebar for structural changes
                var sidebar = document.querySelector('.sidebar-container') || 
                             document.querySelector('aside.sidebar') ||
                             document.querySelector('aside');
                if (sidebar) {
                    sidebarObserver.observe(sidebar, {
                        childList: true,
                        subtree: true
                    });
                }
            }
        }, 300);
    }
    
    // Listen for theme changes
    function setupThemeListener() {
        // Furo theme switcher buttons
        var themeButtons = document.querySelectorAll('button[data-theme-toggle], [data-theme-toggle]');
        themeButtons.forEach(function(button) {
            button.addEventListener('click', function() {
                setTimeout(function() {
                    switcherCreated = false;
                    createLanguageSwitcher();
                }, 300);
            });
        });
        
        // Watch for data-theme attribute changes on body (Furo's way)
        if (!observer) {
            observer = new MutationObserver(function(mutations) {
                var themeChanged = false;
                mutations.forEach(function(mutation) {
                    if (mutation.type === 'attributes' && 
                        (mutation.attributeName === 'data-theme' || mutation.attributeName === 'class')) {
                        themeChanged = true;
                    }
                });
                
                if (themeChanged) {
                    setTimeout(function() {
                        switcherCreated = false;
                        createLanguageSwitcher();
                    }, 200);
                } else {
                    // Check if switcher was removed
                    var switcherExists = document.querySelector('#language-switcher');
                    if (!switcherExists || !switcherExists.parentNode) {
                        switcherCreated = false;
                        setTimeout(createLanguageSwitcher, 100);
                    }
                }
            });
        }
        
        // Observe body for theme changes
        if (document.body) {
            observer.observe(document.body, {
                attributes: true,
                attributeFilter: ['data-theme', 'class']
            });
        }
    }
    
    document.addEventListener('DOMContentLoaded', setupThemeListener);
    
    // Also set up immediately if DOM is already ready
    if (document.readyState !== 'loading') {
        setupThemeListener();
    }
    
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

