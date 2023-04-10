module ConditionalProbability
using DataFrames
using Makie, Makie.GeometryBasics
using Random
using CSV

export ProbProb, CardsPlot, ExportCSV, ExportImages

"""
    ProbProb(A⁺C⁺, A⁺C⁻, A⁻C⁺, A⁻C⁻)
"""
function ProbProb(A⁺C⁺::a, A⁺C⁻::b, A⁻C⁺::c, A⁻C⁻::d) where 
    {a <: Number, b <: Number, c <: Number, d <: Number}
    A⁺ = A⁺C⁺ + A⁺C⁻       #  $A = AC + A\neg C$
    A⁻ = A⁻C⁺ + A⁻C⁻       # $\neg A = \neg AC + \neg A \neg C$
    C⁺ = A⁺C⁺ + A⁻C⁺       #  $C = AC + \neg A C$
    C⁻ = A⁺C⁻ + A⁻C⁻       # $\neg C = A \neg C + \neg A \neg C$
    (Σ = A⁺ + A⁻) == (C⁺ + C⁻) || error("总数不相等")

    C⁺⎸A⁺ = A⁺C⁺ / A⁺       # Conditional Probability: $P(C|A)$
    C⁺⎸A⁻ = A⁻C⁺ / A⁻       # Conditional Probability: $P(C|\neg A)$
    C⁻⎸A⁺ = A⁺C⁻ / A⁺       # Conditional Probability: $P(\neg C|A)$
    MI    = (Σ - A⁺C⁻) / Σ  # Material Implication

    ΔP  = C⁺⎸A⁺ - C⁺⎸A⁻     # $\Delta P = P(C|A) - P(C|\neg A)$
    ΔP⁺ = ΔP / (1 - C⁺⎸A⁻)  # $\Delta P^+ = \frac{\Delta P}{1 - P(C|\neg A)}$ Causal Probality

    dt = [A⁺C⁺ A⁺C⁻ A⁻C⁺ A⁻C⁻ A⁺ C⁺ Σ MI C⁻⎸A⁺ C⁺⎸A⁺ C⁺⎸A⁻ ΔP ΔP⁺]
    df = DataFrame(dt, ["A⁺C⁺", "A⁺C⁻", "A⁻C⁺", "A⁻C⁻", "A⁺", "C⁺", 
    "Σ", "MI", "C⁻|A⁺", "C⁺|A⁺", "C⁺|A⁻", "ΔP", "ΔP⁺"])

    [df[!, x] = convert.(Int, df[!, x])      for x in 1:7]
    [df[!, x] = round.(df[!, x], digits = 2) for x in 8:ncol(df)]

    return df
end

function ProbProb(A, B, C, D)
    dt = [ProbProb(a, b, c, d) for a in A for b in B for c in C for d in D]
    return reduce(vcat, dt)
end

"""
    CardsPlot(A, B, C, D; shuffle = true, sentence = "如果这张卡片上的图形是红色的，那么它是圆形。")
"""
function CardsPlot(A, B, C, D; shuffle = true, 
	sentence = "如果卡片上的图形是红色的，那么它是圆形。"
)
    
    SpaceT = 0.05
    ShiftX = 0.15
    FigW = 6
    FigH = 4
    FigResW = 1000 

    order = shuffle ? randperm(4) .- 1 : [3, 2, 1, 0]

    function Card(shape, color, RecB, n)
        RecL = ShiftX * (n - 1)
        RecH = 1 - SpaceT
        RecW = RecH * 3 / 4

        CircleC = Point2f(RecL + RecW / 2, RecB + RecH / 2)
        CircleR = RecW / 2 - ShiftX

        CubeL = RecL + RecW / 2 - CircleR
        CubeB = RecB + RecH / 2 - CircleR

        poly!(
			ax, Rect(RecL, RecB, RecW, RecH), 
            strokewidth = 2, color = :gray90
		)
        if     shape == :circle
            poly!(ax, Circle(CircleC, CircleR), color = color)
        elseif shape == :square
            poly!(ax, Rect(CubeL, CubeB, 2 * CircleR, 2 * CircleR), 
                color = color)
        end
    end

    fig = Figure(
		resolution = (FigResW, FigResW * FigH / FigW), 
		backgroundcolor = :gray99
	)
	Label(fig[0, 1], sentence, fontsize = 25, 
		padding = (0, 0, -30, -20), tellwidth = false)
    ax = Axis(fig[1, 1], 
        limits = ((-SpaceT, FigW), (-SpaceT, FigH)))
    # hidespines!(ax)
    hidedecorations!(ax)
    [Card(:circle,  :red, order[1], n) for n in 1:A]
    [Card(:square,  :red, order[2], n) for n in 1:B]
    [Card(:circle, :blue, order[3], n) for n in 1:C]
    [Card(:square, :blue, order[4], n) for n in 1:D]

    fig
end

"""
    ExportCSV(data::AbstractDataFrame; dataDir = "Conditional_Materials")
"""
function ExportCSV(data::AbstractDataFrame; dataDir = "Conditional_Materials")
    isdir(dataDir) || mkdir(dataDir)
    DataCSV = transform(data,
        ["A⁺C⁺", "A⁺C⁻", "A⁻C⁺", "A⁻C⁻"] => ByRow(
            (A, B, C, D) -> "P_" * lpad(A, 2, "0") * "_" * lpad(B, 2, "0") * 
            "_" * lpad(C, 2, "0") * "_" * lpad(D, 2, "0") * ".png"
        ) => "Image"
    )
    rename!(DataCSV, Dict("A⁺C⁺" => "A1C1", "A⁺C⁻" => "A1C0", "A⁻C⁺" 
    => "A0C1", "A⁻C⁻" => "A0C0", "A⁺" => "A1", "C⁺" => "C1", "Σ" => 
    "Total", "MI" => "MI", "C⁻|A⁺" => "C0GvA1", "C⁺|A⁺" => "C1GvA1", 
    "C⁺|A⁻" => "C1GvA0", "ΔP" => "dltP", "ΔP⁺" => "dltPp"))
    CSV.write(joinpath(dataDir, "Conditional_Stimuli.csv"), DataCSV)
end

"""
    ExportImages(data::AbstractDataFrame; dataDir = "Conditional_Materials")
"""
function ExportImages(data::AbstractDataFrame; dataDir = "Conditional_Materials")
    isdir(dataDir) || mkdir(dataDir)
    for row in eachrow(data)
        A, B, C, D = row."A⁺C⁺", row."A⁺C⁻", row."A⁻C⁺", row."A⁻C⁻"
        fig = CardsPlot(A, B, C, D)
        save(joinpath(dataDir, "P_" * lpad(A, 2, "0") * "_" * 
        lpad(B, 2, "0") * "_" * lpad(C, 2, "0") * "_" * 
        lpad(D, 2, "0") * ".png"), fig, px_per_unit = 10)
    end
    println(joinpath(pwd(), dataDir))
end

function ExportImages(A::a, B::b, C::c, D::d; dataDir = "Conditional_Materials")  where 
    {a <: Number, b <: Number, c <: Number, d <: Number}
    isdir(dataDir) || mkdir(dataDir)
    fig = CardsPlot(A, B, C, D)
    save(joinpath(dataDir, "P_" * lpad(A, 2, "0") * "_" * 
    lpad(B, 2, "0") * "_" * lpad(C, 2, "0") * "_" * 
    lpad(D, 2, "0") * ".png"), fig, px_per_unit = 10)
    println(joinpath(pwd(), dataDir))
end

end # module
