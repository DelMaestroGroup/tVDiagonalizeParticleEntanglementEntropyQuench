using Base

"""Data type to handle file outputs"""
mutable struct FileOutputHandler
    files_dict::Dict{String, Union{IOStream,Nothing}}
    paths_dict::Dict{String, String} 
    data_to_file_functions_dict::Dict{String, Function} 
    hname_format_dict::Dict{String, String} 
    nEntries::Int64
    flush_on_write::Bool   
end
"""Initialize an empty file output handler"""
FileOutputHandler(flush_on_write::Bool) = FileOutputHandler(Dict{String, IOStream}(),Dict{String, String}(),Dict{String, Function}(),Dict{String, String}(),0,flush_on_write)
FileOutputHandler() = FileOutputHandler(Dict{String, IOStream}(),Dict{String, Function}(),Dict{String, String}(),0,true)
function add!(fh::FileOutputHandler,path::String,data_to_file_function::Function,handler_name::String;hname_format::String="", replace_list::Any=[nothing]) 
    if ~(hname_format == "") 
        fh.hname_format_dict[handler_name] = hname_format
        handler_name = format(handler_name*hname_format,replace_list...)
    end
    for (key,val) in fh.paths_dict
        if key == handler_name
            error("File handeler names must be unique. Tried to add new file handler with already existing name ", handler_name,".")
        end
        if val == path
            error("File with path '", path,"' is already part of the file handler.")
        end
    end 
    fh.nEntries += 1
    fh.paths_dict[handler_name] = path
    fh.data_to_file_functions_dict[handler_name] = data_to_file_function 
    # open the file later only when needed
    fh.files_dict[handler_name] = nothing
    return nothing
end
"""Get file from handler name"""
function get_file(fh::FileOutputHandler,handler_name::String;open_as::String ="w",replace_list::Any=[nothing])
    if haskey(fh.hname_format_dict,handler_name)
        hname_format = fh.hname_format_dict[handler_name]
        handler_name = format(handler_name*hname_format,replace_list...)
    end
    if fh.files_dict[handler_name] === nothing
        fh.files_dict[handler_name] = open(fh.paths_dict[handler_name],open_as)
    end 
    return fh.files_dict[handler_name]
end
"""Write data to file. File is chosen by handler_name, data to string function is the correponding 
function in data_to_file_functions_dict."""
function Base.write(fh::FileOutputHandler,handler_name::String,data...;replace_list::Any=[nothing]) 
    if haskey(fh.hname_format_dict,handler_name)
        hname_format = fh.hname_format_dict[handler_name]
        handler_name = format(handler_name*hname_format,replace_list...)
    end
    if fh.files_dict[handler_name] === nothing 
        fh.files_dict[handler_name] = open(fh.paths_dict[handler_name],"w")
    end 
    write_str = fh.data_to_file_functions_dict[handler_name](data...)
    write_flush(fh.files_dict[handler_name],write_str,fh.flush_on_write)
    return nothing
end
"""Write string to file. File is chosen by handler_name."""
function write_str(fh::FileOutputHandler,handler_name::String,str::String;replace_list::Any=[nothing])  
    if haskey(fh.hname_format_dict,handler_name)
        hname_format = fh.hname_format_dict[handler_name]
        handler_name = format(handler_name*hname_format,replace_list...)
    end
    file = get_file(fh,handler_name;open_as="w")
    write_flush(file,str,fh.flush_on_write)
    return nothing
end
"""Close all files."""
function Base.close(fh::FileOutputHandler) 
    for file in fh.files_dict.vals
        if ~(file === nothing)
            close(file)
        end
    end 
    return nothing
end
function Base.close(fh::FileOutputHandler,handler_name::String;replace_list::Any=[nothing]) 
    if haskey(fh.hname_format_dict,handler_name)
        hname_format = fh.hname_format_dict[handler_name]
        handler_name = format(handler_name*hname_format,replace_list...)
    end
    close(fh.files_dict[handler_name])
    return nothing
end

"""Use 'write' to write string to IOstream (e.g. write to a file) and flush IOstream if toflush is true."""
function write_flush(stream::IO,str::String,toflush::Bool=true)
    write(stream, str)
    if toflush
        flush(stream)
    end
    return nothing
end