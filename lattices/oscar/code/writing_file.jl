LAT_INCLUSIONS_HEADER = (
  ("domain", "text"),
  ("codomain", "text"),
  ("primitive", "boolean"),
  ("quotient_invs", "integer[]"),
  ("domain_rank", "smallint"),
  ("codomain_rank", "smallint"),
  ("quotient_rank", "smallint"),
  ("multiplicity", "integer"),
  ("orthogonal_complement", "numeric[]"),
  ("quotient_gens", "numeric[]"),
  ("images", "numeric[]")
)

default_table_path = joinpath(@__DIR__, "lat_inclusions.tbl")
default_data_path = joinpath(@__DIR__, "../data")
k3_label = "22.16.1.0.1"

function write_data_to_table(
  table_path::String = default_table_path,
  data_path::String = default_data_path,
  header_and_types = LAT_INCLUSIONS_HEADER,
)
  touch(table_path)
  header = join(first.(header_and_types), "|")*"\n"*join(last.(header_and_types), "|")*"\n\n"
  f = open(table_path, "w")
  write(f, header)
  close(f)
  
  for _s in readdir(data_path; join=true)
    i = parse(Int, last(_s))
    for d_label in readdir(_s), _data in readdir(joinpath(_s, d_label); join=true)
      if contains(_data, "k3_embed")
        line_data = [d_label, k3_label, "T", "{}", "$i", "22", "$(22-i)", only(readlines(_data)), "{}", "{}", "{}"]
        line = join(line_data, "|")*"\n\n"
        f = open(table_path, "a")
        write(f, line)
        close(f)
      elseif contains(_data, "overlat")
        _overlats = Dict{NTuple{2, String}, Vector{String}}()
        for _line in readlines(_data)
          c_label, qi, qg = split(_line, "_")
          qi = "{$(qi)}"
          key = (c_label, qi)
          qg = replace(qg, "[" => "{", "]" => "}")
          if haskey(_overlats, key)
            push!(_overlats[key], qg)
          else
            _overlats[key] = String[qg]
          end
        end
        for ((c_label, qi), qgs) in _overlats
          qg = "{"*join(qgs, ",")*"}"
          line_data = [d_label, c_label, "F", qi, "$i", "$i", "0", "$(length(qgs))", "{}", qg, "{}"]
          line = join(line_data, "|")*"\n\n"
          f = open(table_path, "a")
          write(f, line)
          close(f)
        end
      else
        @assert contains(_data, "prim_emb")
        _prim_embs = Dict{String, Vector{String}}()
        for _line in readlines(_data)
          c_label, oc = split(_line, "_")
          oc = replace(oc, "[" => "{", "]" => "}")
          if haskey(_prim_embs, c_label)
            push!(_prim_embs[c_label], oc)
          else
            _prim_embs[c_label] = String[oc]
          end
        end
        for (c_label, ocs) in _prim_embs
          oc = "{"*join(ocs, ",")*"}"
          line_data = [d_label, c_label, "T", "{}", "$i", "$(i+1)", "1", "$(length(ocs))", oc, "{}", "{}"]
          line = join(line_data, "|")*"\n\n"
          f = open(table_path, "a")
          write(f, line)
          close(f)
        end
      end
    end
  end
  return nothing
end
