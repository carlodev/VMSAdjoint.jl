"""
    verifykey(params::Dict{Symbol,Any},keyname; val = false)

It check if the dictionary params has the entry `keyname`. If not it adds the new entry with the value val. It is used to add default values
"""
function verifykey(params::Dict{Symbol,Any},keyname; val = false)
    if !haskey(params, keyname)
        merge!(params,Dict(keyname=>val))
    end
end

function updatekey(params::Dict{Symbol,Any},keyname::Symbol, val)
    if !haskey(params, keyname)
        merge!(params,Dict(keyname=>val))
    else
        params[keyname] = val
    end
end