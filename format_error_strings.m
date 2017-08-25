


%% This function formats avg and std to the same precision keeping extra digits for errors lying in the interval [0.10; 0.13] * 10^x

function [str_mean, str_std] = format_error_strings(flt_mean, flt_std)

% Initialize
str_mean = '';
str_std = '';
int_keep_extra_digit = 0;

% Determine the scale of fluctuations
sgn = sign(flt_std);
expnt = fix(log10(abs(flt_std)));
mant = sgn * 10^(log10(abs(flt_std))-expnt);
if abs(mant) < 1
    mant = mant * 10;
    expnt = expnt - 1;
end
[mant, expnt] = get_mant_and_exp(flt_std);

% Determine if an extra digit should be kept (for 0.10 -- 0.13)
if mant >= 1 && mant <1.35
    int_keep_extra_digit = 1;
end;

% Format the std string 
mant;
int_keep_extra_digit;
num_digits = -expnt + int_keep_extra_digit;
if num_digits >= 0
    format = sprintf('%%.%if', num_digits);
else
    format = '%.0f';
    flt_std = round(flt_std, num_digits);
end;
str_std = sprintf(format, flt_std);

% Format the avg string
% Be sure to keep at least one significant digit
[~, expnt] = get_mant_and_exp(flt_mean);
num_digits = max([- expnt, num_digits]);
flt_mean = round(flt_mean, num_digits);
format = sprintf('%%.%if', num_digits);
str_mean = sprintf(format, flt_mean);


end





%% == A helper function extracting the mantissa and the exponent of a number ==

function [mant, expnt] = get_mant_and_exp(numb)
    sgn = sign(numb);
    expnt = fix(log10(abs(numb)));
    mant = sgn * 10^(log10(abs(numb))-expnt);
    if abs(mant) < 1
        mant = mant * 10;
        expnt = expnt - 1;
    end
end










