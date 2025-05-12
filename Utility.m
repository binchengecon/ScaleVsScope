
classdef Utility
    methods (Static)

        function prob = P_X(p, x)
            prob = (p.alpha / x)^(p.theta); 
        end


        function f_x = f_X(p, x)
            f_x = p.theta * p.alpha^p.theta * x.^( -p.theta - 1); % f_X function
        end

        function f_y = f_Y(p, y)
            f_y = p.theta * y.^( -p.theta - 1); % f_Y function
        end

        function vec = forward_diff(vec)

            vec = [vec(:,2:end) - vec(:,1:end-1); vec(:,end) - vec(:,end-1)];
        end

        function vec = forward_diff_Omega(vec)
            vec = [vec(2:end) - vec(1:end-1); 0];
        end

        function vec = backward_diff(vec)
            vec = [vec(2) - vec(1); vec(2:end) - vec(1:end-1)];
        end

    end

end