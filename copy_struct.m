function new_struct = copy_struct(old_struct)
    for fn = fieldnames(old_struct)'
        new_struct.(fn{1}) = old_struct.(fn{1});
    end
end