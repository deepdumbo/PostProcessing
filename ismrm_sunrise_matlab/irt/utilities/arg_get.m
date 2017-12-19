 function [arg, ii] = arg_get(list, what, default)
%function [arg, ii] = arg_get(list, what, default)
% extract an argument (e.g., nx or pixel_size, etc.)
% from a list (vertically concatenated strings)
% e.g., as generated by arg_pair() or strvcat()
% the *first* match is used.
% Copyright 2005, Jeff Fessler, The University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if nargin == 1 && streq(list, 'test'), arg_get_test, return, end

list(list == '	') = ' ';	% replace tabs with spaces
what = [what ' '];
for ii=1:size(list,1)
	if streq(list(ii,:), what, length(what))
		arg = list(ii,(1+length(what)):end);
		return
	end
end
if isvar('default') && ~isempty(default)
	if isnumeric(default)
		arg = sprintf('%d', default);
	else
		arg = default;
	end
	ii = 0;
else
	fail('cannot find "%s"', what)
end


function arg_get_test
tmp = strvcat('na 9', 'nb 8', 'nxyz 7');
tmp = arg_get(tmp, 'na');
tmp = str2num(tmp);
jf_equal(tmp, 9)
