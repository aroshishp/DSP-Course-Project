function [y_out, y_state, x_state] = apply_cascade_notch(x_n, y_state, x_state, b_all, a_all)
    num_tones = size(b_all, 1);
    x_input = x_n;

    for k = 1:num_tones
        b = b_all(k,:);
        a = a_all(k,:);

        x1 = x_state(k,1);
        x2 = x_state(k,2);

        y1 = y_state(k,1);
        y2 = y_state(k,2);

        y_n = b(1)*x_input + b(2)*x1 + b(3)*x2 - a(2)*y1 - a(3)*y2;

        x_state(k,2) = x_state(k,1);
        x_state(k,1) = x_input;

        y_state(k,2) = y_state(k,1);
        y_state(k,1) = y_n;

        x_input = y_n;
    end

    y_out = x_input;
end
