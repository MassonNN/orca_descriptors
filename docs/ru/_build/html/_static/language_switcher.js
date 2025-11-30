// Language switcher for ORCA Descriptors documentation

(function() {
    'use strict';
    
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
        
        // Replace language in path
        var newPath = currentPath;
        
        if (currentPath.includes('/' + currentLang + '/')) {
            newPath = currentPath.replace('/' + currentLang + '/', '/' + otherLang + '/');
        } else if (currentPath.includes('/' + otherLang + '/')) {
            newPath = currentPath.replace('/' + otherLang + '/', '/' + currentLang + '/');
        } else {
            // If path doesn't contain language, try to determine from structure
            // For Sphinx, we might be in _build/html/ or similar
            var pathParts = currentPath.split('/');
            var langIndex = -1;
            
            // Look for 'en' or 'ru' in path
            for (var i = 0; i < pathParts.length; i++) {
                if (pathParts[i] === 'en' || pathParts[i] === 'ru') {
                    langIndex = i;
                    break;
                }
            }
            
            if (langIndex >= 0) {
                pathParts[langIndex] = otherLang;
                newPath = pathParts.join('/');
            } else {
                // Fallback: add language to path
                // Try to insert after 'docs' or at the beginning
                if (currentPath.startsWith('/docs/')) {
                    newPath = currentPath.replace('/docs/', '/docs/' + otherLang + '/');
                } else {
                    newPath = '/' + otherLang + currentPath;
                }
            }
        }
        
        // Preserve hash if present
        var newUrl = newPath + currentHash;
        window.location.href = newUrl;
    }
    
    function createLanguageSwitcher() {
        // Wait for DOM to be ready
        if (document.readyState === 'loading') {
            document.addEventListener('DOMContentLoaded', createLanguageSwitcher);
            return;
        }
        
        // Wait a bit for theme to load
        setTimeout(function() {
            var currentLang = getCurrentLanguage();
            var otherLang = getOtherLanguage();
            var langNames = {
                'en': 'English',
                'ru': 'Русский'
            };
            
            // Try to find sidebar or header for Furo theme
            var sidebar = document.querySelector('.sidebar-container') || 
                         document.querySelector('aside') ||
                         document.querySelector('.sidebar');
            
            // Create language switcher button
            var switcher = document.createElement('div');
            switcher.id = 'language-switcher';
            
            var button = document.createElement('button');
            button.textContent = langNames[otherLang];
            button.className = 'language-switcher-btn';
            button.onclick = switchLanguage;
            
            switcher.appendChild(button);
            
            // Add to sidebar if found (Furo theme)
            if (sidebar) {
                var sidebarHeader = sidebar.querySelector('.sidebar-brand') || 
                                   sidebar.querySelector('h1') ||
                                   sidebar.firstElementChild;
                if (sidebarHeader && sidebarHeader.nextSibling) {
                    sidebarHeader.parentNode.insertBefore(switcher, sidebarHeader.nextSibling);
                } else {
                    sidebar.insertBefore(switcher, sidebar.firstChild);
                }
            } else {
                // Fallback: add as floating button
                switcher.style.cssText = 'position: fixed; top: 10px; right: 10px; z-index: 1000;';
                document.body.appendChild(switcher);
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
        }, 200);
    }
    
    // Initialize
    createLanguageSwitcher();
})();

