document$.subscribe(({ body }) => { 
  renderMathInElement(body, {
    delimiters: [
      { left: "$$",  right: "$$",  display: false },
      { left: "$",   right: "$",   display: true },
    //   { left: "\\(", right: "\\)", display: true },
    //   { left: "\[", right: "\]", display: true },
      { left: "```math", right: "```", display: false }
    ],
  })
})