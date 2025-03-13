using DataFrames, CSV, Statsbase

cd("light_chain_novel")

caps = readdir()

compil_tab = DataFrame()

for i in caps
    if startswith(i, "CAP")
        inter = readdir(i)
        igk = i*"/"*inter[1]*"/IGK/igdiscover_new/final/expressed_J.tab"
        igl = i*"/"*inter[1]*"/IGL/igdiscover_new/final/expressed_J.tab"
        igk_missing = true
        igl_missing = true
        if isfile(igk)
            igk_file = CSV.read(igk, DataFrame)
            igk_missing = false
        end
        if isfile(igl)
            igl_file = CSV.read(igl, DataFrame)
            igl_missing = false
        end
        if igk_missing | igl_missing
            if igk_missing & igl_missing
                continue
            elseif igk_missing
                light_tab = igl_file
            elseif igl_missing 
                light_tab = igk_file
            end
        else
            light_tab =  vcat(igk_file, igl_file)   
        end
        id = repeat([i], nrow(light_tab))
        light_tab[!, :cap_id] = id
        compil_tab = vcat(compil_tab, light_tab)
    end
end

CSV.write("light_Jallele_expression.csv", compil_tab)

expr_tab = CSV.read("light_Jallele_expression.csv", DataFrame)
igk_tab = filter(:gene => n -> startswith(n, "IGK"), expr_tab)
igl_tab = filter(:gene => n -> startswith(n, "IGL"), expr_tab)
caps = unique(expr_tab.cap_id)
igkgenes = unique(getindex.(split.(unique(igk_tab.gene), "*"), 1))
iglgenes = unique(getindex.(split.(unique(igl_tab.gene), "*"), 1))

igkdict = Dict()
igldict = Dict()



for i in caps
    cap_igk = filter(:cap_id => n -> n == i, igk_tab)
    cap_igl = filter(:cap_id => n -> n == i, igl_tab)
    igk_tot = sum(cap_igk.unique_CDR3)
    igl_tot = sum(cap_igl.unique_CDR3)

    for g in igkgenes
        gene_tab = filter(:gene=> n -> startswith(n, g), cap_igk)
        prop = sum(gene_tab.unique_CDR3)/igk_tot
        if haskey(igkdict, g)
            push!(igkdict[g], prop)
        else
            igkdict[g] = [prop]
        end
    end

    for g in iglgenes
        gene_tab = filter(:gene=> n -> startswith(n, g), cap_igl)
        prop = sum(gene_tab.unique_CDR3)/igl_tot
        if haskey(igldict, g)
            push!(igldict[g], prop)
        else
            igldict[g] = [prop]
        end
    end
end

CSV.write("igkJ_expr_prop.csv", DataFrame(igkdict))
CSV.write("iglJ_expr_prop.csv", DataFrame(igldict))

