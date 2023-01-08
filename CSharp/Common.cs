using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using JA.Dynamics;
using JA.LinearAlgebra;

namespace JA
{
    public static class Common
    {
        public static string VariableName(string symbol, string subscript = null)
        {
            if (!string.IsNullOrEmpty(subscript))
            {
                return $"{symbol}_{subscript}";
            }
            return symbol;
        }
        public static string VariableNames(string symbol, string components, string subscript = null)
        {
            var parts = components.Split(',');
            for (int i = 0; i < parts.Length; i++)
            {
                parts[i] = VariableName($"{symbol}{parts[i].Trim()}", subscript);
            }
            return string.Join(",", parts);
        }

    }
}
