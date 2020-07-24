using Documenter, multiflap

makedocs(doctest = false,
	sitename = "Multiple shootin algorithm documentation",
	format = Documenter.HTML(collapselevel = 1),
	# format = DocumenterLaTeX.LaTeX(),
	authors = "Gianmarco Ducci",
	pages = Any[
		"Home" => "index.md",
		"Tutorials" => "tutorials.md",
	]
	)

deploydocs(
	repo = "https://github.com/vortexlab-uclouvain/multiflap.git",
)
