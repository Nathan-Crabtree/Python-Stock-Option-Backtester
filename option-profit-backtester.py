# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.

import QuantLib as ql
import pandas_datareader as pdr
import datetime as dt

from math import *


def main(ticker, volatility, option_time_period, percent_investment_per_trade, initial_capital, year, month, day):

    #start_date = dt.datetime(year, month, day)
    start_date = dt.datetime(year-1, month, day)
    end_date = dt.datetime.today().date()
    #print("end_date: ", end_date)

    stock_price = pdr.get_data_yahoo(ticker, start_date, end_date)

    #print(stock_price.head())

    #print(stock_price.iloc[1]['Close'])
    #print(stock_price.iloc[4286]['Close'])
    #print(len(stock_price.index))
    #set the Date column to be accessable using ['Date']
    stock_price = stock_price.reset_index()
    current_total_capital = initial_capital

    #252 is movivng clock forward to account for volatility calculation
    #for i in range(252, len(stock_price.index) - option_time_period - 1, option_time_period):
    for i in range(0, len(stock_price.index) - option_time_period - 1, option_time_period):

        print("Date: ", stock_price.iloc[i]['Date'])
        investment_capital = current_total_capital * percent_investment_per_trade * (1/100)
        print("Current Total Capital: ", current_total_capital)
        print("Investment Capital: ", investment_capital)
        print(ticker, " price today: ", stock_price.iloc[i]['Close'])
        print(ticker, " price on the day before expiration: ", stock_price.iloc[i + option_time_period - 1]['Close'])
        print("Option Expiration Date: ", stock_price.iloc[i + option_time_period]['Date'])

        #print(int(stock_price.iloc[i]['Date'].day),' ', int(stock_price.iloc[i]['Date'].month),' ',int(stock_price.iloc[i]['Date'].year))
        date = ql.Date(int(stock_price.iloc[i]['Date'].day), int(stock_price.iloc[i]['Date'].month), int(stock_price.iloc[i]['Date'].year))

        maturity_date = ql.Date(int(stock_price.iloc[i + option_time_period]['Date'].day),
                                int(stock_price.iloc[i + option_time_period]['Date'].month),
                                int(stock_price.iloc[i + option_time_period]['Date'].year))
        strike_price = stock_price.iloc[i]['Close']
        spot_price = stock_price.iloc[i]['Close']
        # volatility = 0.20  # the historical vols or implied vols
        # dividend_rate = 0.0163
        # option_type = ql.Option.Call
        # risk_free_rate = 0.001
        # day_count = ql.Actual365Fixed()
        # calendar = ql.UnitedStates()
        # calculation_date = ql.Date(int(stock_price.iloc[i]['Date'].day),
        #                            int(stock_price.iloc[i]['Date'].month),
        #                            int(stock_price.iloc[i]['Date'].year))
        # ql.Settings.instance().evaluationDate = calculation_date
        #
        # payoff = ql.PlainVanillaPayoff(option_type, strike_price)
        # settlement = calculation_date
        #
        # am_exercise = ql.AmericanExercise(settlement, maturity_date)
        # american_option = ql.VanillaOption(payoff, am_exercise)
        #
        # # Black-Scholes-Merton
        # spot_handle = ql.QuoteHandle(
        #     ql.SimpleQuote(spot_price)
        # )
        # flat_ts = ql.YieldTermStructureHandle(
        #     ql.FlatForward(calculation_date, risk_free_rate, day_count)
        # )
        # dividend_yield = ql.YieldTermStructureHandle(
        #     ql.FlatForward(calculation_date, dividend_rate, day_count)
        # )
        # flat_vol_ts = ql.BlackVolTermStructureHandle(
        #     ql.BlackConstantVol(calculation_date, calendar, volatility, day_count)
        # )
        # bsm_process = ql.BlackScholesMertonProcess(spot_handle,
        #                                            dividend_yield,
        #                                            flat_ts,
        #                                            flat_vol_ts)
        #
        # steps = 200
        # binomial_engine = ql.BinomialVanillaEngine(bsm_process, "crr", steps)
        # american_option.setPricingEngine(binomial_engine)
        # #print(american_option.NPV())
        # value_when_bought = american_option.NPV()

        value_when_bought = calculate_option_price(date, maturity_date, volatility, spot_price, strike_price)
        #value_when_bought = calculate_option_price(date, maturity_date, calculate_volatility(stock_price, i), spot_price, strike_price)
        total_contracts_bought = investment_capital // (value_when_bought * 100)
        print("Total Contracts Bought: ", total_contracts_bought)



        #calculate option value at expiration

        # maturity_date = ql.Date(int(stock_price.iloc[i + option_time_period]['Date'].day), int(stock_price.iloc[i + option_time_period]['Date'].month), int(stock_price.iloc[i + option_time_period]['Date'].year))
        # strike_price = stock_price.iloc[i]['Close']
        # spot_price = stock_price.iloc[i + option_time_period -1]['Close']
        # volatility = 0.20  # the historical vols or implied vols
        # dividend_rate = 0.0163
        # option_type = ql.Option.Call
        # risk_free_rate = 0.001
        # day_count = ql.Actual365Fixed()
        # calendar = ql.UnitedStates()
        # calculation_date = ql.Date(int(stock_price.iloc[i + option_time_period -1]['Date'].day), int(stock_price.iloc[i + option_time_period -1]['Date'].month), int(stock_price.iloc[i + option_time_period -1]['Date'].year))
        # ql.Settings.instance().evaluationDate = calculation_date
        #
        # payoff = ql.PlainVanillaPayoff(option_type, strike_price)
        # settlement = calculation_date
        #
        # am_exercise = ql.AmericanExercise(settlement, maturity_date)
        # american_option = ql.VanillaOption(payoff, am_exercise)
        #
        # # Black-Scholes-Merton
        # spot_handle = ql.QuoteHandle(
        #     ql.SimpleQuote(spot_price)
        # )
        # flat_ts = ql.YieldTermStructureHandle(
        #     ql.FlatForward(calculation_date, risk_free_rate, day_count)
        # )
        # dividend_yield = ql.YieldTermStructureHandle(
        #     ql.FlatForward(calculation_date, dividend_rate, day_count)
        # )
        # flat_vol_ts = ql.BlackVolTermStructureHandle(
        #     ql.BlackConstantVol(calculation_date, calendar, volatility, day_count)
        # )
        # bsm_process = ql.BlackScholesMertonProcess(spot_handle,
        #                                            dividend_yield,
        #                                            flat_ts,
        #                                            flat_vol_ts)
        #
        # steps = 200
        # binomial_engine = ql.BinomialVanillaEngine(bsm_process, "crr", steps)
        # american_option.setPricingEngine(binomial_engine)
        #print(american_option.NPV())
        # value_at_expiration = american_option.NPV()

        calculation_date = ql.Date(int(stock_price.iloc[i + option_time_period -1]['Date'].day), int(stock_price.iloc[i + option_time_period -1]['Date'].month), int(stock_price.iloc[i + option_time_period -1]['Date'].year))
        maturity_date = ql.Date(int(stock_price.iloc[i + option_time_period]['Date'].day), int(stock_price.iloc[i + option_time_period]['Date'].month), int(stock_price.iloc[i + option_time_period]['Date'].year))
        strike_price = stock_price.iloc[i]['Close']
        spot_price_expiration = stock_price.iloc[i + option_time_period -1]['Close']

        value_at_expiration = calculate_option_price(calculation_date, maturity_date, volatility, spot_price_expiration, strike_price)
       #value_at_expiration = calculate_option_price(calculation_date, maturity_date, calculate_volatility(stock_price, i), spot_price_expiration, strike_price)
        profit = (value_at_expiration - value_when_bought) * 100 * total_contracts_bought

        if value_at_expiration < 0.0001:
            value_at_expiration = 0
        print("Value when bought: ", value_when_bought)
        print("Value at expiration: ", value_at_expiration)
        print("Profit: ", profit)
        print("")
        print("")

        current_total_capital = (current_total_capital + profit)

    date = ql.Date(5, 8, 2006)
    print(date)
    print(date.year())
    print(date.month())
    print(date.dayOfMonth())

    # option data
    maturity_date = ql.Date(29, 8, 2006)
    spot_price = 22.79
    strike_price = 23
    volatility = 1.18  # the historical vols or implied vols
    dividend_rate = 0.0163
    option_type = ql.Option.Call

    risk_free_rate = 0.001
    day_count = ql.Actual365Fixed()
    calendar = ql.UnitedStates()

    calculation_date = ql.Date(5, 8, 2006)
    ql.Settings.instance().evaluationDate = calculation_date

    payoff = ql.PlainVanillaPayoff(option_type, strike_price)
    settlement = calculation_date

    am_exercise = ql.AmericanExercise(settlement, maturity_date)
    american_option = ql.VanillaOption(payoff, am_exercise)

    #Black-Scholes-Merton
    spot_handle = ql.QuoteHandle(
        ql.SimpleQuote(spot_price)
    )
    flat_ts = ql.YieldTermStructureHandle(
        ql.FlatForward(calculation_date, risk_free_rate, day_count)
    )
    dividend_yield = ql.YieldTermStructureHandle(
        ql.FlatForward(calculation_date, dividend_rate, day_count)
    )
    flat_vol_ts = ql.BlackVolTermStructureHandle(
        ql.BlackConstantVol(calculation_date, calendar, volatility, day_count)
    )
    bsm_process = ql.BlackScholesMertonProcess(spot_handle,
                                               dividend_yield,
                                               flat_ts,
                                               flat_vol_ts)

    steps = 200
    binomial_engine = ql.BinomialVanillaEngine(bsm_process, "crr", steps)
    american_option.setPricingEngine(binomial_engine)
    print(american_option.NPV())
    print(bs_call(spot_price, strike_price, 0.1, 0.001, 0.2))

def calculate_option_price(calculation_date, maturity_date, volatility, spot_price, strike_price):

    # print("Calculation Date: ", calculation_date)
    # print("Maturity Date: ", maturity_date)
    # print("Spot Price: ", spot_price)
    # print("Strike Price: ", strike_price)

    #volatility = 0.20  # the historical vols or implied vols
    dividend_rate = 0.0163
    option_type = ql.Option.Call

    risk_free_rate = 0.001
    day_count = ql.Actual365Fixed()
    calendar = ql.UnitedStates()

    ql.Settings.instance().evaluationDate = calculation_date

    payoff = ql.PlainVanillaPayoff(option_type, strike_price)
    settlement = calculation_date

    am_exercise = ql.AmericanExercise(settlement, maturity_date)
    american_option = ql.VanillaOption(payoff, am_exercise)

    # Black-Scholes-Merton
    spot_handle = ql.QuoteHandle(
        ql.SimpleQuote(spot_price)
    )
    flat_ts = ql.YieldTermStructureHandle(
        ql.FlatForward(calculation_date, risk_free_rate, day_count)
    )
    dividend_yield = ql.YieldTermStructureHandle(
        ql.FlatForward(calculation_date, dividend_rate, day_count)
    )
    flat_vol_ts = ql.BlackVolTermStructureHandle(
        ql.BlackConstantVol(calculation_date, calendar, volatility, day_count)
    )
    bsm_process = ql.BlackScholesMertonProcess(spot_handle,
                                               dividend_yield,
                                               flat_ts,
                                               flat_vol_ts)

    steps = 200
    binomial_engine = ql.BinomialVanillaEngine(bsm_process, "crr", steps)
    american_option.setPricingEngine(binomial_engine)
    return(american_option.NPV())

def calculate_volatility(stock_price, i):
    # print("i: ", i)
    # print("stock price[i]: ", stock_price.iloc[i]['Close'])
    # print("stock price[i-251]: ", stock_price.iloc[i-251]['Close'])
    # print("stock price[i-252]: ", stock_price.iloc[i-252]['Close'])

    count = 0
    sum_interday_value_change_percentage = 0
    sum_interday_value_change_percentage_minus_mean_squared = 0
    sum_logarithmic_returns = 0
    for j in range(i-251, i):
        #sum_interday_value_change_percentage += (stock_price.iloc[j]['Close'] / stock_price.iloc[j-1]['Close'] - 1)
        sum_logarithmic_returns += log(stock_price.iloc[j]['Close'] / stock_price.iloc[j - 1]['Close'])
        count += 1
    print
    #mean = sum_interday_value_change_percentage / count
    mean = sum_logarithmic_returns / count
    #print("Mean: ", mean)
    count = 0
    sum_logarithmic_returns_minus_mean_squared = 0
    for j in range(i-251, i):
        #sum_interday_value_change_percentage_minus_mean_squared += ((stock_price.iloc[j]['Close'] / stock_price.iloc[j-1]['Close'] - 1) - mean) ** 2
        sum_logarithmic_returns_minus_mean_squared += ((log(stock_price.iloc[j]['Close']) / stock_price.iloc[j - 1]['Close']) - mean) ** 2
        count +=1
    #print("Volatility: ", ((1/count)*sum_interday_value_change_percentage_minus_mean_squared) ** (1/2))
    print("Volatility: ", ((1/count) * sum_logarithmic_returns_minus_mean_squared) ** (1/2))
   # return ((1/count)*sum_interday_value_change_percentage_minus_mean_squared) ** (1/2)
    return ((1/count) * sum_logarithmic_returns_minus_mean_squared) ** (1/2)

#first define these 2 functions
def d1(S,X,T,r,sigma):
    return (log(S/X)+(r+sigma*sigma/2.)*T)/(sigma*sqrt(T))

def d2(S,X,T,r,sigma):
    return d1(S,X,T,r,sigma)-sigma*sqrt(T)

#define the call option price function
def bs_call(S,X,T,r,sigma):
     return S*CND(d1(S,X,T,r,sigma))-X*exp(-r*T)*CND(d2(S,X,T,r,sigma))

#define the put options price function
def bs_put(S,X,T,r,sigma):
      return X*exp(-r*T)-S + bs_call(S,X,T,r,sigma)

#define cumulative standard normal distribution
def CND(X):
     (a1,a2,a3,a4,a5)=(0.31938153,-0.356563782,1.781477937,-1.821255978,1.330274429)
     L = abs(X)
     K=1.0/(1.0+0.2316419*L)
     w=1.0-1.0/sqrt(2*pi)*exp(-L*L/2.)*(a1*K+a2*K*K+a3*pow(K,3)+a4*pow(K,4)+a5*pow(K,5))
     if X<0:
        w=1.0-w
     return w

if __name__ == '__main__':
    main(ticker='AMZN', volatility = 0.3, option_time_period = 60, percent_investment_per_trade = 30, initial_capital = 10000, year = 2006, month = 8, day = 1)


