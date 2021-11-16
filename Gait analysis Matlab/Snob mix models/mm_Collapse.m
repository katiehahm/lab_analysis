function mm = mm_Collapse(mm, ix)

%% if all classes are to be purged...
if(all(ix))
    ix(1) = false;  % save one class 
end

mm.nClasses = mm.nClasses - sum(ix);        % remaining classes

mm.a(ix)    = [];                           % new mixing proportions
mm.a        = mm.a ./ sum(mm.a);            % normalize so that sum a_i = 1
mm.class    = removecells(mm.class, ix);    % remove class details

end

