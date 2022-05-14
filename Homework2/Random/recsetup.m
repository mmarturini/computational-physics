function s = recsetup(a,b)
    s = b;
    for fname = fieldnames(a)'
        if isstruct(a.(fname{1}))
            s.(fname{1}) = recsetup(a.(fname{1}),b.(fname{1}));
        else
            s.(fname{1}) = a.(fname{1});
        end
    end
end