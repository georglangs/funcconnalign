function p = defaultParams(p,defaultp)

names = fieldnames(defaultp);

for n = 1:length(names)
    name = names{n};

    if not(isfield(p,name))
        p = setfield(p,name,getfield(defaultp,name));
    end
end