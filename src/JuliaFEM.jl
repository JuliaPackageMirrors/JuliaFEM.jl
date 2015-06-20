# This file is a part of JuliaFEM. License is MIT: https://github.com/ovainola/JuliaFEM/blob/master/README.md
module JuliaFEM

export Model, new_model, new_field, get_field, add_nodes, get_nodes


@doc """
Basic model
""" ->
type Model
    model  # For global variables
    nodes  # For nodes
    elements  # For elements
    element_nodes  # For element nodes
    element_gauss_points  # For element gauss points
end




@doc """
Initialize empty model.

Parameters
----------
None

Returns
-------
New model struct
""" ->
function new_model()
    return Model(Dict(), Dict(), Dict(), Dict(), Dict())
end




@doc """
Create new field to model.

Parameters
----------
field_type : Dict()
    Target topology (model.model, model.nodes, model.elements,
                     model.element_nodes, model.element_gauss
field_name : str
    Field name

Returns
-------
Dict
    New field

""" ->
function new_field(field_type, field_name)
    d = Dict()
    setindex!(field_type, d, field_name)
    return d
end




@doc """Get field from model.

Parameters
----------
field_type : Dict()
    Target topology (model.model, model.nodes, model.elements,
                     model.element_nodes, model.element_gauss
field_name : str
    Field name

create_if_doesnt_exist : bool, optional
    If field doesn't exists, create one and return empty field

Raises
------
Error, if field not found and create_if_doesnt_exist == false
""" ->
function get_field(field_type, field_name; create_if_doesnt_exist=false)
    if !(field_name in keys(field_type))
        if create_if_doesnt_exist
            field = new_field(field_type, field_name)
        else
            throw("Field not found")
        end
    else
        field = getindex(field_type, field_name)
    end
    return field
end




@doc """Add new nodes to model.

Parameters
----------
nodes : Dict()
    id => coords

Returns
-------
model

Notes
-----
Create new field "coords" to model if not found
""" ->
function add_nodes(model, nodes)
    println("Adding ", length(nodes), " nodes to model")
    field = get_field(model.nodes, "coords"; create_if_doesnt_exist=true)
    for (node_id, coords) in nodes
        field[node_id] = coords
    end
end




@doc """Return subset of nodes from model.

Parameters
----------
node_ids : array
    list of node ids to return

Returns
-------
Dict()
    id => coords
""" ->
function get_nodes(model, node_ids)
    subset = Dict()
    for node_id in node_ids
        subset[node_id] = model.nodes["coords"][node_id]
    end
    return subset
end


end # module