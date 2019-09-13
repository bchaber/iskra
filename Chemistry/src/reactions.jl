import DataStructures: OrderedDict

macro reactions(ex)
  def_reaction_network(ex)
end

struct Reaction
  rate # might be any symbol: function or variable
  reactants # list of reactant's name and its order
  stoichiometry # list of net coefficients 
end

function def_reaction_network(ex::Expr)
  reactions = :([])
  mapping = OrderedDict{String,Int}()

  for arg in ex.args
    parse(arg, reactions, mapping)
  end

  return reactions
end

function parse(ex::LineNumberNode, reactions, mapping) end

function parse(ex::Expr, reactions, mapping)
    if ex.head == :tuple
      rate = ex.args[1]
      eqn  = ex.args[2]

      reactants = :([])     # indices
      stoichiometry = :([]) # net stoichiometry

      reacs, prods = parse_reaction(eqn, mapping)
      for (id, ix) in mapping
        P = get(prods, string(id), 0)
        R = get(reacs, string(id), 0)
        if haskey(reacs, id)
          push!(reactants.args, :((string($id), $R)))
        end

        net = P - R

        if net != 0
          push!(stoichiometry.args, :((string($id), $net)))
        end
      end

      push!(reactions.args, :(Reaction($(esc(rate)), $reactants, $stoichiometry)))
    end
end

function parse_reaction(ex::Expr, mapping)
  reactants = Dict{String,Int}()
  products  = Dict{String,Int}()

  if ex.head == :-->
    exr = ex.args[1] # LHS
    exp = ex.args[2] # RHS

    add_participants!(reactants, exr, mapping)
    add_participants!(products,  exp, mapping)
  else
    throw("malformed reaction")
  end

  return reactants, products
end

function add_participants!(dict, ex, mapping)
  if isa(ex, Symbol) # found a species symbol
    if !haskey(mapping, ex)
      mapping[string(ex)] = length(mapping) + 1
    end

    id = ex
    val = get(dict, string(id), 0)

    dict[string(id)] = val + 1
  
  elseif isa(ex, Expr) && ex.args[1] == :* # species symbol has a coefficient
    id    = ex.args[3]
    coeff = ex.args[2]

    if !haskey(mapping, id)
      mapping[string(id)] = length(mapping) + 1
    end

    val = get(dict, string(id), 0)
    dict[string(id)] = val + coeff

  elseif isa(ex, Expr) # found something else, probably of the form (a A + b B)
    for i in 2:length(ex.args)
      add_participants!(dict, ex.args[i], mapping)
    end
  end
end