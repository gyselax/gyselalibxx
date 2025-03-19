window.MathJax = {
  tex: {
    inlineMath: [["$", "$"],["\\(", "\\)"], ["$`", "`$"]],
    displayMath: [["\\[", "\\]"], ["```math", "```"]],
    // processEscapes: true,
    // processEnvironments: true
  }
};

document$.subscribe(() => { 
  MathJax.startup.output.clearCache()
  MathJax.typesetClear()
  MathJax.texReset()
  MathJax.typesetPromise()
})


window.MathJax = {
    startup: {
      ready: () => {
        console.log('MathJax is loaded, but not yet initialized');
        MathJax.startup.defaultReady();
        console.log('MathJax is initialized, and the initial typeset is queued');
      }
    }
  };