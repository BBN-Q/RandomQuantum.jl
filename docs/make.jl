using Documenter, RandomQuantum

makedocs(
	sitename = "RandomQuantum",
	modules  = [RandomQuantum],
	format   = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"),
	pages    = [
		"index.md"
	]
)

