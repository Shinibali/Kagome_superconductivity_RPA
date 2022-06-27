function [bareIntmat] = full_bare_int(nsite,norb,U,J)

Up = U; Jp = J;
[~,~,~,~,bare_int_ind] = basis_formation(nsite,norb); dimen = numel(bare_int_ind);

rows=zeros(1,1) ; cols=rows ; vals=rows ;
for isp2 = 1:2
    for isp1 = 1:2
        for s1 = 1:nsite
            for o2 = 1:norb
                for o1 = 1:norb
                    if o1 == o2 && isp2 ~= isp1
                        value = U;
                        
                        rw = bare_int_ind( o1,s1,o1,s1,isp2,isp1 );
                        cl = bare_int_ind( o1,s1,o1,s1,isp1,isp2 );
                        rows = vertcat(rows, rw(:) ); % joint 1&2
                        cols = vertcat(cols, cl(:) ); % joint 3&4
                        vals = vertcat(vals, value .* ones(size(rw(:))) );
                        
                        rw = bare_int_ind( o1,s1,o1,s1,isp2,isp2 );
                        cl = bare_int_ind( o1,s1,o1,s1,isp1,isp1 );
                        rows = vertcat(rows, rw(:) ); % joint 1&2
                        cols = vertcat(cols, cl(:) ); % joint 3&4
                        vals = vertcat(vals, -value .* ones(size(rw(:))) );
                        
                    elseif o1 ~= o2 && isp2 == isp1
                        value = Up-J ;
                        
                        rw = bare_int_ind( o1,s1,o2,s1,isp2,isp1 );
                        cl = bare_int_ind( o2,s1,o1,s1,isp1,isp2 );
                        rows = vertcat(rows, rw(:) ); % joint 1&2
                        cols = vertcat(cols, cl(:) ); % joint 3&4
                        vals = vertcat(vals, value .* ones(size(rw(:))) );
                        
                        rw = bare_int_ind( o1,s1,o1,s1,isp2,isp2 );
                        cl = bare_int_ind( o2,s1,o2,s1,isp1,isp1 );
                        rows = vertcat(rows, rw(:) ); % joint 1&2
                        cols = vertcat(cols, cl(:) ); % joint 3&4
                        vals = vertcat(vals, -value .* ones(size(rw(:))) );
                        
                    elseif o1 ~= o2 && isp2 ~= isp1
                        value = Up ;
                        
                        rw = bare_int_ind( o1,s1,o2,s1,isp2,isp1 );
                        cl = bare_int_ind( o2,s1,o1,s1,isp1,isp2 );
                        rows = vertcat(rows, rw(:) ); % joint 1&2
                        cols = vertcat(cols, cl(:) ); % joint 3&4
                        vals = vertcat(vals, value .* ones(size(rw(:))) );
                        
                        rw = bare_int_ind( o1,s1,o1,s1,isp2,isp2 );
                        cl = bare_int_ind( o2,s1,o2,s1,isp1,isp1 );
                        rows = vertcat(rows, rw(:) ); % joint 1&2
                        cols = vertcat(cols, cl(:) ); % joint 3&4
                        vals = vertcat(vals, -value .* ones(size(rw(:))) );
                        
                        value = J ;
                        
                        rw = bare_int_ind( o1,s1,o1,s1,isp2,isp1 );
                        cl = bare_int_ind( o2,s1,o2,s1,isp1,isp2 );
                        rows = vertcat(rows, rw(:) ); % joint 1&2
                        cols = vertcat(cols, cl(:) ); % joint 3&4
                        vals = vertcat(vals, value .* ones(size(rw(:))) );
                        
                        rw = bare_int_ind( o1,s1,o2,s1,isp2,isp2 );
                        cl = bare_int_ind( o2,s1,o1,s1,isp1,isp1 );
                        rows = vertcat(rows, rw(:) ); % joint 1&2
                        cols = vertcat(cols, cl(:) ); % joint 3&4
                        vals = vertcat(vals, -value .* ones(size(rw(:))) );
                        
                        value = Jp ;
                        
                        rw = bare_int_ind( o1,s1,o2,s1,isp2,isp1 );
                        cl = bare_int_ind( o1,s1,o2,s1,isp1,isp2 );
                        rows = vertcat(rows, rw(:) ); % joint 1&2
                        cols = vertcat(cols, cl(:) ); % joint 3&4
                        vals = vertcat(vals, value .* ones(size(rw(:))) );
                        
                        rw = bare_int_ind( o1,s1,o2,s1,isp2,isp2 );
                        cl = bare_int_ind( o1,s1,o2,s1,isp1,isp1 );
                        rows = vertcat(rows, rw(:) ); % joint 1&2
                        cols = vertcat(cols, cl(:) ); % joint 3&4
                        vals = vertcat(vals, -value .* ones(size(rw(:))) );
                        
                    end
                end
            end
        end
    end
end
rows=rows(2:end) ; cols=cols(2:end) ; vals=vals(2:end) ;
bareIntmat = sparse(rows,cols,vals,dimen,dimen);

end