## Script to analyse UK and NOR sequences ##

# A. Marsden - 2024 #

################################################################################################################

#defining functions

using CSV, DataFrames, ProgressBars, HypothesisTests, StatsBase, Statistics

function compile_haplo(directory)

    p_dirs = readdir(directory*"/P1")

    p_dirs_filtered = []
    for i in p_dirs
        if split(i, "_")[1] == "P1"
            push!(p_dirs_filtered,i)
        end
    end

    bar_dirs = ProgressBar(p_dirs_filtered)

    gene_list = []

    for i in p_dirs_filtered
        p_geno = CSV.read(open(directory*"/P1/"*i*"/"*i*"_genotype.tsv"), DataFrame, delim="\t")
        for x in p_geno.gene
            if x ∉ gene_list
                push!(gene_list, x)
            end
        end
    end



    part_dict = Dict()

    

    for i in bar_dirs
        print(i)
        alleles = []
        for gene in gene_list
            p_geno = CSV.Rows(open(directory*"/P1/"*i*"/"*i*"_genotype.tsv"), delim="\t")
            all_missing = true
            for r in p_geno
                if r.gene == gene
                    push!(alleles, r.GENOTYPED_ALLELES)
                    all_missing = false
                end
            end

            if all_missing
                push!(alleles, missing)
            end
        end

        part_dict[i] = alleles 
    end

    output_tab = DataFrame(part_dict)
    output_tab[!, :gene] = gene_list

    return output_tab

end

function expression_lc(directory)
    p_dirs = readdir(directory*"/P1")

    p_dirs_filtered = []
    for i in p_dirs
        if split(i, "_")[1] == "P1"
            push!(p_dirs_filtered,i)
        end
    end

    gene_list = []

    for i in p_dirs_filtered
        p_geno = CSV.read(open(directory*"/P1/"*i*"/"*i*"_genotype.tsv"), DataFrame, delim="\t")
        for x in p_geno.gene
            if x ∉ gene_list
                if occursin("J",x)
                    continue
                end
                push!(gene_list, x)
            end
        end
    end

    expr_dict = Dict()
    bar_dirs = ProgressBar(p_dirs_filtered)
    for i in bar_dirs
        print(i)
        exprss = []
        p_tab = CSV.read(open(directory*"/P1/"*i*"/"*i*"_genotype.tsv"), DataFrame, delim="\t")
        p_tab = filter(:total => n -> n != "NA", p_tab)
        p_tab = filter(:gene => n -> occursin("J",n), p_tab)
        if typeof(p_tab.total) == Vector{Int64}
            total = sum( p_tab.total)
        else
            total = sum(parse.(Int64, p_tab.total))
        end

        for gene in gene_list
            p_geno = CSV.Rows(open(directory*"/P1/"*i*"/"*i*"_genotype.tsv"), delim="\t")
            all_missing = true
            for r in p_geno
                if r.total == "NA"
                    continue
                end
                if r.gene == gene
                    if typeof(r.total) == Int64
                        push!(exprss, r.total/total)
                    else
                        push!(exprss, parse(Int64,r.total)/total)
                    end
                    all_missing = false
                end
            end

            if all_missing
                push!(exprss, 0)
            end
        end

        expr_dict[i] = exprss

    end

    output_tab = DataFrame(expr_dict)
    output_tab[!, :gene] = gene_list

    return output_tab
    
end



dir_igk = "/Users/alaine/Data/Light_Chain/norwegian_lc/IGK"
dir_igl = "/Users/alaine/Data/Light_Chain/norwegian_lc/IGL"

igk_tab = compile_haplo(dir_igk)
igl_tab = compile_haplo(dir_igl)

## need to redo function with all genes first then iterate through gene_list and make new table

igk_expr = expression_lc(dir_igk)
igl_expr = expression_lc(dir_igl)

CSV.write("/Users/alaine/Data/Light_Chain/norwegian_lc/igk_expr_nor.csv", igk_expr)
CSV.write("/Users/alaine/Data/Light_Chain/norwegian_lc/igl_expr_nor.csv", igl_expr)

nor_igk_expr = CSV.read("/Users/alaine/Data/Light_Chain/norwegian_lc/igk_expr_nor.csv", DataFrame)
nor_igl_expr = CSV.read("/Users/alaine/Data/Light_Chain/norwegian_lc/igl_expr_nor.csv", DataFrame)

cap_igk_expr = CSV.read("/Users/alaine/Data/Light_Chain/re_analysis/igk_expr_prop.csv", DataFrame)
cap_igl_expr = CSV.read("/Users/alaine/Data/Light_Chain/re_analysis/igl_expr_prop.csv", DataFrame)

igk_gene_list = names(cap_igk_expr)
igl_gene_list = names(cap_igl_expr)

test_dict_igk = Dict()

for i in igk_gene_list
    nor_par = 0
    for x in eachrow(nor_igk_expr)
        if i == x.gene
            print("success")
            for x in eachrow(nor_igk_expr)
                if i == x.gene
                    nor_par = collect(x)
                    pop!(nor_par)
                end
            end
        end
    end

    cap_par = cap_igk_expr[!,i]

    runs = []
    for y in 1:1000
        nor_sample = sample(nor_par,6, replace=false)
        test = UnequalVarianceTTest(cap_par, convert(Vector{Float64}, nor_sample))
        push!(runs, pvalue(test))
    end

    test_dict_igk[i] = mean(runs)

end

test_dict_igl = Dict()

for i in igl_gene_list
    nor_par = 0
    for x in eachrow(nor_igl_expr)
        if i == x.gene
            print("success")
            for x in eachrow(nor_igl_expr)
                if i == x.gene
                    nor_par = collect(x)
                    pop!(nor_par)
                end
            end
        end
    end

    cap_par = cap_igl_expr[!,i]

    runs = []
    for y in 1:1000
        nor_sample = sample(nor_par,6, replace=false)
        test = UnequalVarianceTTest(cap_par, convert(Vector{Float64}, nor_sample))
        push!(runs, pvalue(test))
    end

    test_dict_igl[i] = mean(runs)

end

CSV.write("/Users/alaine/Data/Light_Chain/igk_expr_comp.csv", DataFrame(test_dict_igk))
CSV.write("/Users/alaine/Data/Light_Chain/igl_expr_comp.csv", DataFrame(test_dict_igl))


filter!(row -> row[:gene] in igk_gene_list, igk_expr)
filter!(row -> row[:gene] in igl_gene_list, igl_expr)

CSV.write("/Users/alaine/Data/Light_Chain/norwegian_lc/igk_expr_nor.csv", igk_expr)
CSV.write("/Users/alaine/Data/Light_Chain/norwegian_lc/igl_expr_nor.csv", igl_expr)

#### uk repertoire stat generation

rename_dict = Dict(
    "IGKV1D-12" => "IGKV1E-12",
    "IGKV1-12" => "IGKV1E-12",
    "IGKV1D-13" => "IGKV1E-13",
    "IGKV1-13" => "IGKV1E-13",
    "IGKV1D-33" => "IGKV1E-33",
    "IGKV1-33" => "IGKV1E-33",
    "IGKV1D-37" => "IGKV1E-37",
    "IGKV1-37" => "IGKV1E-37",
    "IGKV1D-39" => "IGKV1E-39",
    "IGKV1-39" => "IGKV1E-39",
    "IGKV2D-28" => "IGKV2E-28",
    "IGKV2-28" => "IGKV2E-28",
    "IGKV2D-40" => "IGKV2E-40",
    "IGKV2-40" => "IGKV2E-40",
    "IGKV6D-21" => "IGKV6E-21",
    "IGKV6-21" => "IGKV2E-21",
)

brits = readdir("/Users/alaine/Downloads/airr_files")

igk_expr_dict = Dict()
igl_expr_dict = Dict()

for file in brits
    tab = CSV.read("/Users/alaine/Downloads/airr_files/"*file, DataFrame)
    gene_list = unique(first.(split.(tab.v_call, "*")))
    igk_list = []
    igl_list = []
    for i in gene_list
        if startswith(i, "IGK")
            push!(igk_list, i)
        elseif startswith(i, "IGL")
            push!(igl_list, i)
        end
    end
    igk_tot = 0
    igl_tot = 0
    for row in eachrow(tab)
        if startswith(row.v_call, "IGK")
            igk_tot=igk_tot+1
        elseif startswith(row.v_call, "IGL")
            igl_tot=igl_tot+1
        end
    end


    for i in igk_list
        gene_count = 0
        for x in tab.v_call
            if startswith(x, i)
                gene_count = gene_count+1
            end
        end
        prop = gene_count/igk_tot
        if i ∉ keys(igk_expr_dict)
            igk_expr_dict[i] = [prop]
        else
            push!(igk_expr_dict[i], prop)
        end
    end

    for i in igl_list
        gene_count = 0
        for x in tab.v_call
            if startswith(x, i)
                gene_count = gene_count+1
            end
        end
        prop = gene_count/igl_tot
        if i ∉ keys(igl_expr_dict)
            igl_expr_dict[i] = [prop]
        else
            push!(igl_expr_dict[i], prop)
        end
    end
end


function dict_fill(dict)

    max_val = maximum(length.(values(dict)))

    for (k,v) in dict
        if length(v) == max_val
            continue
        else
            for i in 1:(max_val-length(v))
                push!(v, 0)
            end
        end
    end

end

uk_igk = DataFrame(igk_expr_dict)
uk_igl = DataFrame(igl_expr_dict)

select!(uk_igk, igk_gene_list)
select!(uk_igl, igl_gene_list)

CSV.write("/Users/alaine/Data/Light_Chain/igk_expr_uk.csv", uk_igk)
CSV.write("/Users/alaine/Data/Light_Chain/igl_expr_uk.csv", uk_igl)

test_dict_igk = Dict()

for i in igk_gene_list
    uk_par = uk_igk[!,i]


    cap_par = cap_igk_expr[!,i]

    runs = []
    for y in 1:100
        nor_sample = sample(uk_par,6, replace=false)
        test = UnequalVarianceTTest(cap_par, uk_par)
        push!(runs, pvalue(test))
    end

    test_dict_igk[i] = mean(runs)

end

test_dict_igl = Dict()

for i in igl_gene_list
    uk_par = uk_igl[!,i]


    cap_par = cap_igl_expr[!,i]

    runs = []
    for y in 1:100
        nor_sample = sample(uk_par,6, replace=false)
        test = UnequalVarianceTTest(cap_par, uk_par)
        push!(runs, pvalue(test))
    end

    test_dict_igl[i] = mean(runs)

end

CSV.write("/Users/alaine/Data/Light_Chain/igk_expr_ukxsa_comp.csv", DataFrame(test_dict_igk))
CSV.write("/Users/alaine/Data/Light_Chain/igl_expr_ukxsa_comp.csv", DataFrame(test_dict_igl))

test_dict_igk = Dict()

for i in igk_gene_list
    nor_par = 0
    for x in eachrow(nor_igk_expr)
        if i == x.gene
            print("success")
            for x in eachrow(nor_igk_expr)
                if i == x.gene
                    nor_par = collect(x)
                    pop!(nor_par)
                end
            end
        end
    end

    uk_par = uk_igk[!,i]

    runs = []
    for y in 1:1000
        nor_sample = sample(nor_par,10, replace=false)
        test = UnequalVarianceTTest(uk_par, convert(Vector{Float64}, nor_sample))
        push!(runs, pvalue(test))
    end

    test_dict_igk[i] = mean(runs)

end

test_dict_igl = Dict()

for i in igl_gene_list
    nor_par = 0
    for x in eachrow(nor_igl_expr)
        if i == x.gene
            print("success")
            for x in eachrow(nor_igl_expr)
                if i == x.gene
                    nor_par = collect(x)
                    pop!(nor_par)
                end
            end
        end
    end

    uk_par = uk_igl[!,i]

    runs = []
    for y in 1:1000
        nor_sample = sample(nor_par,10, replace=false)
        test = UnequalVarianceTTest(uk_par, convert(Vector{Float64}, nor_sample))
        push!(runs, pvalue(test))
    end

    test_dict_igl[i] = mean(runs)

end

CSV.write("/Users/alaine/Data/Light_Chain/igk_expr_ukxnor_comp.csv", DataFrame(test_dict_igk))
CSV.write("/Users/alaine/Data/Light_Chain/igl_expr_ukxnor_comp.csv", DataFrame(test_dict_igl))
