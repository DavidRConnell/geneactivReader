function values = readGeneactivBin(fileName)
	fileStr = fileread(fileName);
	[values.startTime, values.sampFreq, values.header, ...
		dataStartIndex, firstTempIndex, calibrationData] = readHeader(fileStr);

	[acc, temp] = readAndConvertBinData(fileStr, dataStartIndex, firstTempIndex);
	values.temp = temp;
	values.acc = calibrateXyz(acc, calibrationData);

	function [startTime, sampFreq, header, ...
		dataStartIndex, firstTempIndex, calibrationData] = readHeader(fileStr)

		startTime = getStartTime(fileStr);
		sampFreq = getSamplingFreq(fileStr);
		[header, dataStartIndex, firstTempIndex] = getHeader(fileStr);
		calibrationData = getCalibrationData(fileStr);

		function startTime = getStartTime(fileStr)
			timeIndex = regexp(fileStr, 'Page Time', 'once');
			startTime = fileStr(timeIndex + 21:timeIndex + 32);
		end

		function sampFreq = getSamplingFreq(fileStr)
			searchPattern = 'Measurement Frequency:[0-9]* Hz';
			freqString = regexp(fileStr, searchPattern, 'match', 'once');
			sampFreq = regexp(freqString, '[0-9]*', 'match', 'once');
			sampFreq = str2double(sampFreq);
		end

		function [header, dataStartIndex, firstTempIndex] = getHeader(fileStr)
			dataStartIndex = regexp(fileStr, 'Recorded Data', 'once');
			header = fileStr(1:(dataStartIndex - 3));
			firstTempIndex = regexp(fileStr, 'Temperature:', 'once', 'end');
		end

		function calibrationData = getCalibrationData(fileStr)
			searchPattern = 'Calibration Data[A-Za-z0-9\n\r:-\s]*Volts';
			calibrationData = regexp(fileStr, searchPattern, 'match', 'once');
			calibrationData = regexp(calibrationData, '[0-9-]*', 'match');
		end
	end
end

function [acc, temp] = readAndConvertBinData(fileStr, dataStartIndex, firstTempIndex)
	findPatternIndex = @(pattern) regexp(fileStr, pattern, 'once');
	charsPerPageInitial = findPatternIndex('Sequence Number:1') - ...
		findPatternIndex('Sequence Number:0');

	pagesStr = regexp(fileStr(dataStartIndex:-1:1), '[0-9]+', 'once', 'match');
	numPages = str2double(pagesStr(end:-1:1));

	if ~iscompletebin(numPages, fileStr)
		error('MATLAB:corruptFile', 'Number of pages not equal to check.')
	end

	pagesOfDataPerLoop = 10;
	samplesPerPage = 300;
	hexCharsPerPage = samplesPerPage * 12;

	startIdx = 1;
	hexStart = dataStartIndex + charsPerPageInitial - 2 - hexCharsPerPage;

	acc = zeros(3, numPages * 300);
	indeces = getIndecesOfAccelDataInHex(pagesOfDataPerLoop);
	hexChars = zeros(1, hexCharsPerPage * pagesOfDataPerLoop);
	temp = zeros(1, numPages);
	currentTempIndex = firstTempIndex + 1;

	for pageNum = 1:pagesOfDataPerLoop:numPages
		if pageNum + pagesOfDataPerLoop - 1 > numPages
			pagesOfDataPerLoop = numPages - pageNum;
			indeces = getIndecesOfAccelDataInHex(pagesOfDataPerLoop);
		end

		for pageStepper = 0:(pagesOfDataPerLoop - 1)
			exactPageNum = pageNum + pageStepper;
			charsPerPage = charsPerPageInitial + floor(log10(exactPageNum));

			[temp(exactPageNum), currentTempIndex, charsPerPage] = ...
				getTemp(currentTempIndex, charsPerPage);

			hexStoreStart = pageStepper * hexCharsPerPage + 1;
			hexChars(hexStoreStart:(hexStoreStart + hexCharsPerPage - 1)) = ...
				fileStr(hexStart:(hexStart + hexCharsPerPage - 1));

			hexStart = hexStart + charsPerPage;
		end

		acc(:, startIdx:(startIdx + 300 * pagesOfDataPerLoop - 1)) = ...
			hex2xyz(hexChars, indeces);

		startIdx = startIdx + 300 * pagesOfDataPerLoop;
	end
	acc = convertMsbToSignBit(acc);

	function flag = iscompletebin(numPages, fileStr)
		pattern = strcat('\d*', reverse('Sequence Number:'));
		lastPageNum = reverse(regexp(reverse(fileStr), pattern, 'once', 'match'));
		lastPageNum = str2double(regexp(lastPageNum, '\d*', 'once', 'match')) + 1;
		flag = numPages == lastPageNum;
	end

	function indeces = getIndecesOfAccelDataInHex(pagesOfDataPerLoop)
		hexCharsPerDatum = 3;
		shapeIndeces = [hexCharsPerDatum, ...
			hexCharsPerPage * pagesOfDataPerLoop / hexCharsPerDatum];

		indeces = reshape(1:(hexCharsPerPage * pagesOfDataPerLoop), shapeIndeces)';

		numSignals = 4;
		numRows = samplesPerPage * numSignals;

		desiredRows = mod(1:(numRows * pagesOfDataPerLoop), 4) ~= 0;
		desiredRows = desiredRows' .* (1:(numRows * pagesOfDataPerLoop))';
		indeces = indeces(desiredRows ~= 0, :);
	end

	function [temp, nextTempIndex, charsPerPage] = ...
		getTemp(currentTempIndex, charsPerPage)

		tempStr = fileStr(currentTempIndex:(currentTempIndex + 3));
		temp = str2double(tempStr);
		charsPerPage = fixCharsPerPage(charsPerPage);
		nextTempIndex = currentTempIndex + charsPerPage;

		function charsPerPage = fixCharsPerPage(charsPerPage)
			expectedStrLength = 4;
			actualStrLength = length(num2str(temp, '%0.1f'));
			charsPerPage = charsPerPage - (expectedStrLength - actualStrLength);
		end
	end

	function acc = hex2xyz(hexString, indeces)
		acc = hex2dec(hexString(indeces));
		acc = reshape(acc, [3 length(acc) / 3]);
	end

	function acc = convertMsbToSignBit(acc)
		N = 12;
		acc = 2 ^ (N - 1) - acc;
		acc = sign(acc) .* (2 ^ (N - 1) - abs(acc));
	end
end

function acc = calibrateXyz(acc, calibrationData)
	[x_gain, x_offset, y_gain, y_offset, z_gain, z_offset] = calibrationData{:};

	acc = acc * 100 - [str2double(x_offset); str2double(y_offset); str2double(z_offset)];
	acc = acc ./ [str2double(x_gain); str2double(y_gain); str2double(z_gain)];
end
