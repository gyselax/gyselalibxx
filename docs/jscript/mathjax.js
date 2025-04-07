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
        console.log('MathJax is loaded, but not yet initialised');
        MathJax.startup.defaultReady();
        console.log('MathJax is initialised, and the initial typeset is queued');
      }
    }
  };