using IJulia
using JSON
using Pkg
ENV["GKSwstype"] = "100"

function main(;
    basedir=get(ENV, "DOCDIR", "docs"),
    cachedir=get(ENV, "NBCACHE", ".cache"),
    rmsvg=true)
    nb = get(ENV, "NB", "test.ipynb")
    IJulia.installkernel("Julia", "--project=@.")
    # nbconvert command options
    kernelname = "--ExecutePreprocessor.kernel_name=julia-1.$(VERSION.minor)"
    execute = ifelse(get(ENV, "ALLOWERRORS", "false") == "true", "--execute --allow-errors", "--execute")
    timeout = "--ExecutePreprocessor.timeout=" * get(ENV, "TIMEOUT", "-1")
    nbout = joinpath(abspath(pwd()), cachedir, nb)
    mkpath(dirname(nbout))
    cmd = `jupyter nbconvert --to notebook $(execute) $(timeout) $(kernelname) --output $(nbout) $(nb)`
    run(cmd)
    rmsvg && strip_svg(nbout)
    return nothing
end

# Strip SVG output from a Jupyter notebook
function strip_svg(ipynb)
    oldfilesize = filesize(ipynb)
    nb = open(JSON.parse, ipynb, "r")
    for cell in nb["cells"]
        !haskey(cell, "outputs") && continue
        for output in cell["outputs"]
            !haskey(output, "data") && continue
            datadict = output["data"]
            if haskey(datadict, "image/png") || haskey(datadict, "image/jpeg")
                delete!(datadict, "text/html")
                delete!(datadict, "image/svg+xml")
            end
        end
    end
    rm(ipynb; force=true)
    write(ipynb, JSON.json(nb, 1))
    @info "Stripped SVG in $(ipynb). The original size is $(Base.format_bytes(oldfilesize)). The new size is $(Base.format_bytes(filesize(ipynb)))."
    return ipynb
end

main()
