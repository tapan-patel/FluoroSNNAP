function h = PleaseWait(msg)
h = figure('units','pixels','position',[500 500 400 80],'windowstyle','modal');
uicontrol('style','text','string',msg,'units','pixels','position',[2 20 400 40]);