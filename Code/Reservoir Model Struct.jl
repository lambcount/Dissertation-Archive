mutable struct ReservoirModel{T<:Number,N} <: AbstractArray{T,N}
    Population::Array{T,N}
    Bleach::Array{T,N}
    Parameter::Vector{T}
    Time
    PumpedMode
    Control
    Dictionary
end

Base.size(model::ReservoirModel) = size(model.Bleach)
Base.getindex(model::ReservoirModel{T,N}, I::Vararg{Int, N}) where {N,T} = getindex(model.Bleach, I...)
Base.setindex!(model::ReservoirModel{T,N}, v::Number, I::Vararg{Int, N}) where {N,T} = setindex!(model.Bleach, v, I...)
Base.copy(model::ReservoirModel) = ReservoirModel(copy(model.Population), copy(model.Bleach),copy(model.Parameter), copy(model.Time),copy(model.PumpedMode),copy(model.Control),copy(model.Dictionary))