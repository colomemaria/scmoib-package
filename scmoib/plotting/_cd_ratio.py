def cd_ratio(df):
    tdf = df.copy()
    tdf['Connected'] = 1 - tdf['disc_ratio']
    tdf.rename(columns={'disc_ratio': 'Disconnected'}, inplace=True)
    tdf[['Disconnected', 'Connected']].plot(kind='bar', stacked=True)
