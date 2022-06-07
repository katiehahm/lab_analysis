function [final_sequences] = dfs_IDsequence2(n, N, seg, impact_times, last_o_time, last_x_time, consecutive_o, consecutive_x, final_sequences)
% 4/20/22 DFS search for possible sequences

    % if end of the 
    if n == N
%         seg
        final_sequences = [final_sequences;seg];
    else

        thru_o = false;
        % add o
        time_diff = impact_times(n) - last_o_time;
        if time_diff < 1.6 && time_diff > 0.1
            thru_o = true;
            seg(end+1) = 1;
            if n > 1
                if last_o_time == impact_times(n-1)
                    consecutive_o = consecutive_o + 1;
                else
                    consecutive_o = 1;
                end
            end
            if consecutive_o < 4
                final_sequences = dfs_IDsequence(n+1,N,seg,impact_times,impact_times(n),last_x_time, consecutive_o, consecutive_x,final_sequences);
            end
        end

        % add x
        time_diff = impact_times(n) - last_x_time;
        if time_diff < 1.6 && time_diff > 0.1
            if thru_o
                seg(end) = 2;
            else
                seg(end+1) = 2;
            end
            if n > 1
                if last_x_time == impact_times(n-1)
                    consecutive_x = consecutive_x + 1;
                else
                    consecutive_x = 1;
                end
                
            end
            if consecutive_x < 4
                final_sequences = dfs_IDsequence(n+1,N,seg,impact_times,last_o_time,impact_times(n),consecutive_o, consecutive_x,final_sequences);
            end
        end
    end
end

