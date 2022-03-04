import DataStructures: OrderedDict

macro reactions(ex)
  def_reaction_network(ex)
end

struct ChemicalReaction
  type # reactiono type
  rate # might be any symbol: function or variable
  reactants # list of reactant's name and its order
  stoichiometry # list of net coefficients
end

function def_reaction_network(ex::Expr)
  reactions = :(ChemicalReaction[])
  mapping = OrderedDict{Symbol,Int}()

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
      type = nothing
      if length(ex.args) > 2
        type = ex.args[3]
      end
      reactants = :([])     # indices
      stoichiometry = :([]) # net stoichiometry

      reacs, prods = parse_reaction(eqn, mapping)
      for (id, ix) in mapping
        P = get(prods, (id), 0)
        R = get(reacs, (id), 0)
        if haskey(reacs, id)
          push!(reactants.args, :(($(esc(id)), $R)))
        end

        net = P - R

        if net != 0
          push!(stoichiometry.args, :(($(esc(id)), $net)))
        end
      end

      push!(reactions.args, :(ChemicalReaction($(esc(type)), $(esc(rate)), $reactants, $stoichiometry)))
    end
end

function parse_reaction(ex::Expr, mapping)
  reactants = Dict{Symbol,Int}()
  products  = Dict{Symbol,Int}()

  if ex.head == :--> || ex.head == :->
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
      mapping[(ex)] = length(mapping) + 1
    end

    id = ex
    val = get(dict, (id), 0)

    dict[(id)] = val + 1
  
  elseif isa(ex, Expr) && ex.args[1] == :* # species symbol has a coefficient
    id    = ex.args[3]
    coeff = ex.args[2]

    if !haskey(mapping, id)
      mapping[(id)] = length(mapping) + 1
    end

    val = get(dict, (id), 0)
    dict[(id)] = val + coeff

  elseif isa(ex, Expr) # found something else, probably of the form (a A + b B)
    for i in 2:length(ex.args)
      add_participants!(dict, ex.args[i], mapping)
    end
  end
end